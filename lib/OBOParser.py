#!/usr/bin/env python
'''
    QTLSearch — to search for candidate causal genes in QTL studies
     by combining Gene Ontology annotations across many species, leveraging
    hierarchical orthologous groups.

    (C) 2015-2018 Alex Warwick Vesztrocy <alex@warwickvesztrocy.co.uk>

    This file is part of QTLSearch. It contains a module for parsing an OBO
    file.

    QTLSearch is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    QTLSearch is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with QTLSearch.  If not, see <http://www.gnu.org/licenses/>.
'''
from collections import defaultdict
from itertools import chain
from functools import lru_cache
from numbers import Number
from property_manager import lazy_property
import numpy as np
import re

VERSIONS_SUPPORTED = frozenset({1.2, 1.4})

# Pattern specifications - move these into a config file.
STANZA_TYPE_PATTERN = re.compile(r'\[(?P<type>[^]]*)\]')
TAG_LINE_PATTERN = re.compile(r'^(?P<tag>.+?): *(?P<value>.+?) '
                              r'?(?P<trailing_modifiers>(?<!\\)\{.*?(?<!\\)\})?'
                              r' ?(?P<comment>(?<!\\)!.*?)?$')

# Singularity specification
SINGULAR_TAGS = {
    'header': {
        'format-version',
        'data-version',
        'version',  # deprecated
        'ontology',
        'date',
        'saved-by',
        'auto-generated-by',
        'default-relationship-id-prefix'},
    'instance': {
        'id',
        'is_anonymous',
        'name',
        'namespace',
        'comment',
        'instance_of',
        'created_by',
        'creation_date',
        'is_obsolete'},
    'typedef': {
        'id',
        'is_anonymous',
        'name',
        'namespace',
        'def',
        'domain',
        'range',
        'is_cyclic',
        'is_reflexive',
        'is_symmetric',
        'is_anti_symmetric',
        'is_transitive',
        'is_metadata_tag',
        'is_class_level',
        'is_obsolete'},
    'term': {
        'id',
        'is_anonymous',
        'name',
        'namespace',
        'def',
        'comment',
        'is_obsolete',
        'builtin',
        # Additional tags in 1.4:
        'created_by',
        'creation_date'}}

# Relationship definitions.
# Note: changed for our use here in UBERON.s
UP_RELS = frozenset(['is_a', 'located_in', 'part_of'])
DOWN_RELS = frozenset(['can_be', 'has_member', 'has_part'])
INV_RELS = {'is_a': 'can_be',
            'can_be': 'is_a',
            'part_of': 'has_part',
            'has_part': 'part_of',
            'located_in': 'has_member',
            'has_member': 'located_in'}


ONTOLOGIES = {}


# Custom exceptions
class ObsoleteOntologyTerm(Exception):
    pass


class UnknownOntologyTerm(Exception):
    pass


class OBO(object):
    '''
        This is an extension of the MultiDiGraph networkx class that parses and
        loads an OBO file into it. Supports a fairly strict understanding of the
        v1.2 and v1.4 of the OBO standard, for single-filed ontology
        definitions.

        Term IDs are optionally stored as integers, with an exception if the
        prefix is not that which was expected.
    '''
    def __init__(self, fn, store_as_int=False, ontology_name=None):
        '''
            Loads the OBO file into memory, so that it can be traversed.
        '''
        # Initialise our OBO object.
        self.header = OBOHeader(self)
        self.stanza = OBOStanza(self)

        # Set up option to store as integer
        self._store_as_int = store_as_int
        if store_as_int:
            self.id_prefix = ''
            self.id_length = 0

        # Read the OBO in.
        self._read_obo(fn)

        # Set the name / date - None required for date for some GOSlims.
        self.name = self.header['ontology']
        self.date = self.header['date'] if 'date' in self.header else None

        # Add self to list of ontologies
        if ontology_name is not None:
            ONTOLOGIES[ontology_name] = self

    def __contains__(self, id):
        '''
            Returns boolean dependentant on the term with id being in the
            ontology.
        '''
        try:
            self.get(id)
            return True
        except UnknownOntologyTerm:
            return False
        except ObsoleteOntologyTerm:
            return False

    def __getitem__(self, id):
        '''
            Returns a term from the ontology, converting to int if necessary
        '''
        return self.get(id)

    def __iter__(self):
        yield from self.stanza._stanza['term'].values()

    def get(self, id):
        '''
            Returns a term from the ontology, converting to int if necessary
        '''
        id_ = self._id_to_int(id)
        if id_ in self.stanza['term']:
            return self.stanza['term'][id_]
        elif id_ in self.stanza.obsolete_terms:
            raise ObsoleteOntologyTerm('Term with ID {} is obsolete.'
                                       .format(id))
        elif id_ in self.stanza.alt_ids:
            return self.stanza['term'][self.stanza.alt_ids[id_]]
        elif id_ in self.stanza.xref_ids:
            return self.stanza['term'][self.stanza.xref_ids[id_]]
        else:
            raise UnknownOntologyTerm('Can\'t find term with ID {} in '
                                      'ontology.'.format(id))

    def _load_id(self, id):
        '''
            Loads ID bits if required for store_as_int. Converts and checks also
            in this case. Else, just returns the ID.
        '''
        if self._store_as_int:
            try:
                (prefix, id_new) = id.split(':')
            except ValueError:
                try:
                    return int(id)
                except ValueError:
                    raise ValueError('Can not convert {} to an integer.'
                                     .format(id))

            if self.id_prefix == '':
                self.id_prefix = prefix
            elif self.id_prefix != prefix:
                raise ValueError('Multiple ID prefixes are not supported when '
                                 'storing as ints.')

            if self.id_length == 0:
                self.id_length = len(id_new)
                self.id_length_same = True
            elif self.id_length != len(id_new):
                self.id_length_same = False

            # Convert to integer
            try:
                id = int(id_new)
            except ValueError:
                raise ValueError('Can not convert {} to an integer.'.format(id))

        return id

    def _id_to_int(self, id):
        '''
            Convert ID to integer if required.
        '''
        if isinstance(id, Number) or not self._store_as_int:
            return id

        else:
            try:
                (_, id) = id.split(':')
            except ValueError:
                pass
            try:
                return int(id)
            except ValueError:
                raise ValueError('Can not convert {} to an integer.'.format(id))

    def _id_to_str(self, id):
        if not isinstance(id, str):
            return '{PREFIX:s}:{ID:0{LENGTH:d}d}'.format(
                PREFIX=self.id_prefix,
                LENGTH=(self.id_length if self.id_length_same else
                        len(str(id))),
                ID=id)
        else:
            return id

    def _ids_to_str(self, ids):
        '''
            Returns a list of ids converted back to strings, if necessary.
        '''
        return [self._id_to_str(id) for id in ids]

    def _read_obo(self, fn):
        '''
            Starts the read in of the OBO file.
        '''
        if type(fn) is str:
            with open(fn, 'rt') as fp:
                self._parse_obo(fp)
        else:
            self._parse_obo(fn)

        self._check_consistency()

        # Add the inverse relationships to each of the terms in the ontology.
        self._add_inverse_rels()

    def _parse_obo(self, fp):
        '''
            Parses OBO, yielding terms.
        '''
        header = True
        for line in fp:
            if line[0] != '!' and line != '\n':
                if not header:
                    # Add stanza
                    self.stanza.add(fp, line)
                else:
                    # Add tag-value to header
                    self.header.add(line)

            elif line == '\n':
                # Completed header
                header = False

    def _check_consistency(self):
        '''
            Removes relationships to terms not defined.
        '''
        defined_terms = set(self.stanza['term'].keys())
        for term in self.stanza['term'].keys():
            self.stanza['term'][term].is_a &= defined_terms
            for rel in self.stanza['term'][term].rels:
                self.stanza['term'][term].rels[rel] &= defined_terms

    def _add_inverse_rels(self):
        '''
            Add the inverse relationships to each term.
        '''
        # Create inverse relationships to add to the other terms
        for t in self.stanza['term'].values():
            for (rel, inv_rel) in INV_RELS.items():
                if rel == 'is_a' or rel in t.rels:
                    for term in (t.rels[rel] if rel != 'is_a' else t.is_a):
                        if inv_rel != 'is_a':
                            self.stanza['term'][term].rels[inv_rel].add(t.id)
                        else:
                            # TODO: fix this so that is_a go into the
                            # relationships dict also to make it more clean.
                            self.stanza['term'][term].is_a.add(t.id)

        # Freeze the terms
        for t in self.stanza['term'].values():
            t.is_a = frozenset(t.is_a)
            for rel in t.rels:
                t.rels[rel] = frozenset(t.rels[rel])

    def parents(self, id, as_str=None, include_self=None):
        '''
            Retrieves all parents of the term recursively, optionally including
            itself.
        '''
        include_self = include_self if include_self is not None else False
        as_str = as_str if as_str is not None else False
        if not isinstance(id, (list, set, frozenset, np.ndarray)):
            pass
        else:
            return frozenset(chain.from_iterable(
                self.parents(i, as_str, include_self) for i in id))

        id_ = self[self._id_to_int(id)].id
        parents = self._parentsR(id_)
        additional = set() if not include_self else {id_}

        if not as_str:
            return frozenset(parents | additional)
        else:
            return frozenset(self._ids_to_str(parents | additional))

    @lru_cache(None)
    def _parentsR(self, id):
        '''
            Retrieves all parents of the term, non-recursive
            Note: ensures that loops aren't encountered.
        '''
        parents = set()
        todo = [p for p in self._parents(id)]

        while len(todo) > 0:
            parent = todo.pop()
            parents.add(parent)
            todo += [p
                     for p in self._parents(parent)
                     if p not in parents]

        return parents

    def _parents(self, id):
        '''
            Retrieves all parents of the term.
        '''
        term = self.get(id)
        # Get UP_RELS
        parents = {t
                   for rel in (UP_RELS & set(term.rels.keys()))
                   for t in term.rels[rel]}

        # Add is_a terms
        parents.update(term.is_a)

        return frozenset(parents)

    def children(self, id, as_str=None, include_self=None):
        '''
            Retrieves all children of the term recursively, optionally including
            itself.
        '''
        if not isinstance(id, (list, set, frozenset, np.ndarray)):
            pass
        else:
            return frozenset(chain.from_iterable(
                self.children(i, as_str, include_self) for i in id))

        id_ = self[self._id_to_int(id)].id
        children = self._childrenR(id_)
        additional = set() if not include_self else {id_}

        if not as_str:
            return frozenset(children | additional)
        else:
            return frozenset(self._ids_to_str(children | additional))

    @lru_cache(None)
    def _childrenR(self, id):
        '''
            Retrieves all children of the term, non-recursive
            so as python-friendly.......
            Note: ensures that loops aren't encountered.
        '''
        children = set()
        todo = [c for c in self._children(id)]

        while len(todo) > 0:
            child = todo.pop()
            children.add(child)
            todo += [c
                     for c in self._children(child)
                     if c not in children]

        return children

    def _children(self, id):
        '''
            Retrieves all children of the term.
        '''
        term = self.get(id)
        # Get DOWN_RELS
        return frozenset({t
                          for rel in (DOWN_RELS & set(term.rels.keys()))
                          for t in term.rels[rel]})

    def best_common_term(self, terms, negatives=None):
        '''
            Retrieve the most recent common term of a set of terms.
        '''
        function = self.parents if not negatives else self.children
        # Get initial candidates
        terms = list(terms)
        candidates = set(function(terms[0], include_self=True))
        # Update candidates with intersection of parents for each term.
        for term in terms[1:]:
            candidates.intersection_update(function(term, include_self=True))

        # Order candidates (common terms)
        if candidates:
            return sorted(candidates,
                          key=lambda x: self.get(x).min_depth,
                          reverse=(not negatives))[0]
        else:
            return None


class OBOLine(object):
    '''
        Represents a tag-value line.
    '''
    def __init__(self, line):
        self._match = TAG_LINE_PATTERN.match(line.rstrip())
        if self._match:
            pass
        else:
            raise ValueError('Expected tag-value pair line.')

    @property
    def tag(self):
        '''
            Gets tag as lower case.
        '''
        return self._match.group('tag').lower()

    @property
    def value(self):
        '''
            Gets the value of a tag-value line.
        '''
        return self._match.group('value')

    @property
    def modifiers(self):
        '''
            Gets trailing modifiers, stripping unwanted characters.
        '''
        return self._match.group('trailing_modifiers').strip('{}')

    @property
    def comments(self):
        '''
            Gets comments of term.
        '''
        return self._match.group('comment').lstrip('! ')


class OBOStanza(object):
    '''
        OBO stanza definition.
    '''
    def __init__(self, OBO):
        self._OBO = OBO
        self._stanza = defaultdict(dict)
        self.obsolete_terms = set()
        self.alt_ids = dict()
        self.xref_ids = dict()

    def __getitem__(self, stz_type):
        return self._stanza[stz_type]

    def add(self, fp, line):
        '''
            Adds a stanza to the stanza object.
        '''
        # Get stanza type
        match = STANZA_TYPE_PATTERN.match(line)
        if match:
            stz_type = match.group('type').lower()
            try:
                id = OBOLine(next(fp)).value
                is_obsolete = False
                alt_ids = []
                stz = {}
                while True:
                    line = next(fp)
                    if line != '\n':
                        e = OBOLine(line)
                        if (stz_type != 'term') or (e.tag != 'is_obsolete' and
                                                    e.tag != 'alt_id'):
                            if e.tag in SINGULAR_TAGS[stz_type]:
                                stz[e.tag] = e.value
                            else:
                                stz.setdefault(e.tag, []).append(e.value)
                        elif (e.tag == 'is_obsolete' and
                                e.value.lower() == 'true'):
                            is_obsolete = True
                            while line != '\n':
                                line = next(fp)
                            break
                        elif e.tag == 'alt_id':
                            alt_ids.append(e.value)
                    else:
                        break
            except StopIteration:
                pass
                # raise ValueError('Unexpected EOF.')

            if stz_type == 'term':
                if not is_obsolete:
                    t = OBOTerm(self._OBO, stz, id)
                    self._stanza[stz_type][t.id] = t

                    # Alternative IDs
                    for id_ in alt_ids:
                        self.alt_ids[self._OBO._load_id(id_)] = t.id

                    # Map external ids
                    if 'xref' in stz:
                        for id_ in stz['xref']:
                            self.xref_ids[self._OBO._load_id(id_)] = t.id

                else:
                    self.obsolete_terms.add(id)
            else:
                self._stanza[stz_type][id] = stz

        else:
            raise ValueError('Expected stanza type line.')


class OBOHeader(object):
    '''
        OBO header definition.
    '''
    def __init__(self, OBO):
        self._OBO = OBO
        self._header = {}

    def __contains__(self, k):
        return k in self._header

    def __getitem__(self, k):
        return self._header[k]

    def add(self, line):
        '''
            Adds a tag-value line to the header.
        '''
        e = OBOLine(line)
        tag = e.tag
        if tag in SINGULAR_TAGS['header']:
            self._header[tag] = e.value
            if tag != 'format-version':
                pass
            elif float(e.value) not in VERSIONS_SUPPORTED:
                raise ValueError('OBO format {:f} is not supported.'
                                 .format(e.value))
        else:
            self._header.setdefault(e.tag, []).append(e.value)


class OBOTerm(object):
    '''
        Defines a single term in the ontology.
    '''
    def __init__(self, OBO, stz, id):
        self._OBO = OBO
        self.id = self._OBO._load_id(id)
        self.name = stz.get('name')
        self.definition = stz.get('def')

        # We should parse this a bit more...
        self.is_a = {OBO._id_to_int(x) for x in stz.pop('is_a', [])}
        self._load_rels(stz)

    def __str__(self):
        '''
            Return ID as a string.
        '''
        return self._OBO._id_to_str(self.id)

    def __int__(self):
        '''
            Return ID as an integer.
        '''
        return self._OBO._id_to_int(self.id)

    def _load_rels(self, stz):
        '''
            Load the relationships into a structure and return them.
        '''
        # Load relationships into dictionary
        self.rels = defaultdict(set)
        for rel in stz.pop('relationship', []):
            (typedef, target_term) = rel.split(' ')
            self.rels[typedef].add(self._OBO._id_to_int(target_term))

    @property
    def parentsR(self):
        '''
            Get parents recursively.
        '''
        return self._OBO.parents(self.id)

    @property
    def parents(self):
        '''
            Get immediate parents.
        '''
        return self._OBO._parents(self.id)

    @property
    def childrenR(self):
        '''
            Get children recursively.
        '''
        return self._OBO.children(self.id)

    @property
    def children(self):
        '''
            Get immediate children.
        '''
        return self._OBO._children(self.id)

    @lazy_property
    def min_depth(self):
        '''
            Get the minimum depth of term.
        '''
        # See if we have any parents. If not, level=0
        parents = self.parents
        if parents:
            # recurse on parents
            depths = set()
            for parent in parents:
                depths.add(self._OBO.get(parent).min_depth)

            # Level is the minimum distance to the root
            return 1 + min(depths)
        else:
            return 0

    @lazy_property
    def max_depth(self):
        '''
            Get the maximum depth of term.
        '''
        # See if we have any parents. If not, depth=0
        parents = self.parents
        if parents:
            # recurse on parents
            depths = set()
            for parent in parents:
                depths.add(self._OBO.get(parent).max_depth)

            # Level is the minimum distance to the root
            return 1 + max(depths)
        else:
            return 0

    @lazy_property
    def aspect(self):
        '''
            Gets the aspect (i.e., Biological Process, Molecular Function or
            Cellular Component) of the term, assuming that we are dealing with
            the GO.
        '''
        parents = self._OBO.parents(self.id, include_self=True)
        if parents & {8150, 'GO:0008150'}:
            # Biological Process
            return 'P'
        elif parents & {3674, 'GO:0003674'}:
            # Molecular Function
            return 'F'
        else:
            # Cellular Component
            return 'C'
