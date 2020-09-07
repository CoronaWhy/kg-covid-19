import rdflib
from rdflib import Namespace
import os
import gzip
import shutil

# Base directory of the kg-covid-19 repository
basedir = "./"
param_dict = {
    "GO": {"nodes_id_raw": 0, "nodes_name_raw": 1, "nodes_iri_raw": 4, "nodes_description_raw": 8,
           "nodes_altlab_raw": 7,
           "edges_s_raw": 0, "edges_o_raw": 2, "edges_p_raw": 3,
           "edges_tsv": basedir + "data/transformed/ontologies/go-plus_edges.tsv",
           "nodes_tsv": basedir + "data/transformed/ontologies/go-plus_nodes.tsv"
           },
    "HP": {"nodes_id_raw": 0, "nodes_name_raw": 1, "nodes_iri_raw": 4, "nodes_description_raw": 6,
           "nodes_altlab_raw": 7,
           "edges_s_raw": 0, "edges_o_raw": 2, "edges_p_raw": 3,
           "edges_tsv": basedir + "data/transformed/ontologies/hp_edges.tsv",
           "nodes_tsv": basedir + "data/transformed/ontologies/hp_nodes.tsv"
           }
}
destdir = basedir + "data/rdf/"
savefreq = 10000

skosns = rdflib.namespace.SKOS
dcns = rdflib.namespace.DC
cord19docsns = Namespace("https://cord19.documents/")
rdfns = rdflib.namespace.RDF

# The tsv's have abbreviated URIs
expansionDictionary = {
    "MESH": Namespace("http://id.nlm.nih.gov/mesh/2016/"),
    "GO": Namespace("http://purl.obolibrary.org/obo/"),
    "SO": Namespace("http://purl.obolibrary.org/obo/"),
    "PO": Namespace("http://purl.obolibrary.org/obo/"),
    "CORD": cord19docsns,
    "OBO": Namespace('http://purl.obolibrary.org/obo/'),
    "CHEBI": Namespace("http://purl.obolibrary.org/obo/"),
    "HP": Namespace("http://purl.obolibrary.org/obo/"),
    "NCBITaxon": Namespace("http://purl.obolibrary.org/obo/"),
    "PR": Namespace("http://purl.obolibrary.org/obo/"),
    "UBERON": Namespace("http://purl.obolibrary.org/obo/"),
    "CL": Namespace("http://purl.obolibrary.org/obo/"),
    "BFO": Namespace("http://purl.obolibrary.org/obo/"),
    "FAO": Namespace("http://purl.obolibrary.org/obo/")
}

OBO_spaces = set([k for k, v in expansionDictionary.items() if v == Namespace("http://purl.obolibrary.org/obo/")])


def makeuri(curie):
    splt = curie.split(":")
    prefix = splt[0]
    sufix = splt[1]
    sufix = sufix.replace("body_", "paragraph_")
    uri = None
    if prefix in OBO_spaces:
        a = curie.replace(":", "_")
        uri = expansionDictionary["OBO"][a]
    else:
        uri = expansionDictionary[prefix][sufix]

    return uri


for pk, pv in param_dict.items():
    print("\n-------------------------\n",pk,"Starting")
    source_edges_tsv = pv["edges_tsv"]
    source_nodes_tsv = pv["nodes_tsv"]

    known_docs = set()
    known_paras = set()
    known_entities = set()
    known_prefixes = set()

    concept_labels = dict()

    # Nodes ----------------------------
    with open(source_nodes_tsv) as fin:
        for rown, row in enumerate(fin):
            if rown < 1:
                continue
            rows = row.strip().split("\t")
            if len(rows) < 3:
                continue
            id_raw = rows[pv["nodes_id_raw"]]
            name_raw = rows[pv["nodes_name_raw"]]
            iri_raw = rows[pv["nodes_iri_raw"]]
            description_raw = rows[pv["nodes_description_raw"]] if len(rows) >= pv["nodes_description_raw"] + 1 else ""
            altlab_raw = rows[pv["nodes_altlab_raw"]] if len(rows) >= pv["nodes_altlab_raw"] + 1 else ""
            uri = rdflib.URIRef(iri_raw)
            if uri is None:
                continue
            desc = rdflib.Literal(description_raw, lang="en") if len(description_raw) > 0 else None
            if len(altlab_raw) > 0:
                alt = [rdflib.Literal(xx, lang="en") for xx in altlab_raw.split("|")]
            else:
                alt = []
            concept_labels[uri] = (rdflib.Literal(name_raw, lang="en"), alt, desc)

    destg = rdflib.Graph()

    graphn = 0
    addedtriples = 0
    concepts_to_add = set()
    # EDGES ----------------------------
    with open(source_edges_tsv) as fin:
        for rown, row in enumerate(fin):
            if rown < 1:
                continue
            rows = row.strip().split("\t")
            s_raw = rows[pv["edges_s_raw"]]
            o_raw = rows[pv["edges_o_raw"]]
            p_raw = rows[pv["edges_p_raw"]]

            if p_raw != "rdfs:subClassOf":
                #print(pk,"\tskipping unknown predicate", p_raw)
                continue
            try:
                small_uri = rdflib.URIRef(makeuri(s_raw))
                large_uri = rdflib.URIRef(makeuri(o_raw))
            except:
                continue
            addedtriples += 1
            destg.add((small_uri, skosns['broader'], large_uri))
            destg.add((large_uri, skosns['narrower'], small_uri))
            concepts_to_add.add(large_uri)
            concepts_to_add.add(small_uri)

            if addedtriples % savefreq == savefreq - 1:
                fn = os.path.join(destdir, pk + "_" + str(graphn) + ".ttl")
                destg.serialize(fn, format="ttl")
                graphn += 1
                # destg = rdflib.Graph()

    fn = os.path.join(destdir, pk + ".ttl")
    print(pk,"\tWhole edges graph had  ", len(destg), " triples")
    destg.serialize(fn, format="ttl")
    graphn += 1

    alllabels = set(concept_labels.keys())
    toadd = alllabels & concepts_to_add
    print(pk,"\thas in total",len(toadd),"concepts")
    for uri in toadd:
        destg.add((uri, rdfns["type"], skosns['Concept']))
        preflab, altlabl, desc = concept_labels[uri]
        destg.add((uri, skosns['prefLabel'], preflab))
        for altlab in altlabl:
            destg.add((uri, skosns['altLabel'], altlab))
        if desc is not None:
            destg.add((uri, skosns["definition"], desc))

    fn = os.path.join(destdir, pk + ".ttl")
    destg.serialize(fn, format="ttl")
    print(pk,"\tAfter adding concept info we have ", len(destg), "triples")
    graphn += 1

    with open(fn, 'rb') as f_in:
        with gzip.open(fn.replace("rdf","skos")+".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)