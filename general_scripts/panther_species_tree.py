#!/usr/bin/env python
'''
    Script to download the PANTHER species tree from their API service.
'''
from dendropy import Tree, Node, Taxon
import json
import urllib


def get_panther_species_tree():
    def generate_tree(z, p=None):
        n = Node(label=z['name'])
        # add the information as in nhx

        if 'children' in z:
            # internal
            n = Node()
            n.annotations['S'] = z['name']
            if 'taxon_id' in z:
                n.annotations['TaxonID'] = z['taxon_id']
            for c in z['children']['node']:
                n.add_child(generate_tree(c, p=n))
        else:
            n = Node(taxon=Taxon(z['short_name']))

        return n
    url = 'http://pantherdb.org/services/oai/pantherdb/speciestree'
    req = urllib.request.urlopen(url)
    z = json.loads(req.read().decode('ascii'))['species_tree']
    return Tree(seed_node=generate_tree(z['node']))

t = get_panther_species_tree()
print(t.as_string('newick', annotations_as_nhx=True, suppress_annotations=False))
