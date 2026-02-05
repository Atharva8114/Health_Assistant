from datetime import datetime

def build_payload(
    text,
    pmid,
    title,
    journal,
    year,
    authors,
    section,
    chunk_index,
    api_query,
    entities=None  
):
    # If no entities were found/passed, use the empty template
    if entities is None:
        entities = {
            "drugs": [],
            "conditions": [],
            "biomarkers": [],
            "symptoms": []
        }

    payload = {
        "text": text,
        "pmid": pmid,
        "title": title,
        "journal": journal,
        "year": year,
        "authors": authors,
        "section": section,
        "chunk_index": chunk_index,
        "source": "pubmed_api",
        "api_query": api_query,
        "retrieved_at": datetime.utcnow().isoformat() + "Z",

        # Now we insert the actual extracted entities
        "entities": entities,

        "relations": [],

        "kg_node_ids": {
            "drugs": [],
            "conditions": [],
            "biomarkers": [],
            "symptoms": []
        },

        "study_type": None,
        "confidence_level": None
    }

    return payload



# from datetime import datetime


# def build_payload(
#     text,
#     pmid,
#     title,
#     journal,
#     year,
#     authors,
#     section,
#     chunk_index,
#     api_query
# ):
#     payload = {
#         "text": text,
#         "pmid": pmid,
#         "title": title,
#         "journal": journal,
#         "year": year,
#         "authors": authors,
#         "section": section,
#         "chunk_index": chunk_index,
#         "source": "pubmed_api",
#         "api_query": api_query,
#         "retrieved_at": datetime.utcnow().isoformat() + "Z",

#         # placeholders – you can later plug NER / KG here
#         "entities": {
#             "drugs": [],
#             "conditions": [],
#             "biomarkers": [],
#             "symptoms": []
#         },

#         "relations": [],

#         "kg_node_ids": {
#             "drugs": [],
#             "conditions": [],
#             "biomarkers": [],
#             "symptoms": []
#         },

#         "study_type": None,
#         "confidence_level": None
#     }

#     return payload
