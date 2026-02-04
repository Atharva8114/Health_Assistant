import uuid
from qdrant_client.models import PointStruct

from chunker import simple_chunk
from embedding import embed_texts
from schema_builder import build_payload
from qdrant_store import QdrantClient

from qdrant_store import get_client, create_collection_if_not_exists, create_indexes, COLLECTION
# rename file if needed
# if you name file qdrant_client.py, import carefully to avoid conflict


def ingest_document(
    document_text,
    pmid,
    title,
    journal,
    year,
    authors,
    section,
    api_query
):
    client = get_client()

    create_collection_if_not_exists(client)
    create_indexes(client)

    chunks = simple_chunk(document_text)

    vectors = embed_texts(chunks)

    points = []

    for i, (chunk, vector) in enumerate(zip(chunks, vectors)):
        payload = build_payload(
            text=chunk,
            pmid=pmid,
            title=title,
            journal=journal,
            year=year,
            authors=authors,
            section=section,
            chunk_index=i,
            api_query=api_query
        )

        point_id = str(uuid.uuid4())


        points.append(
            PointStruct(
                id=point_id,
                vector=vector,
                payload=payload
            )
        )

    client.upsert(
        collection_name=COLLECTION,
        points=points
    )

    print(f"Inserted {len(points)} chunks.")


if __name__ == "__main__":

    document = """
    In patients with type 2 diabetes and chronic kidney disease, metformin use
    was associated with improved glycaemic control but increased risk of lactic
    acidosis in advanced renal impairment. This study evaluated long-term
    safety outcomes in a multicenter cohort.
    """

    ingest_document(
        document_text=document,
        pmid="38723421",
        title="Safety and effectiveness of metformin in patients with chronic kidney disease",
        journal="Journal of Clinical Endocrinology",
        year=2024,
        authors=["A. Smith", "B. Kumar", "C. Lee"],
        section="Abstract",
        api_query="metformin chronic kidney disease heart rate variability"
    )
