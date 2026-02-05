import uuid

from qdrant_client.models import PointStruct

from pubmed_fetcher import fetch_pubmed_abstracts
from chunker import simple_chunk
from embedding import embed_texts
from schema_builder import build_payload
from qdrant_store import (
    get_client,
    create_collection_if_not_exists,
    create_indexes,
    COLLECTION
)


def ingest_from_pubmed(query, max_results=5):

    client = get_client()

    create_collection_if_not_exists(client)
    create_indexes(client)

    papers = fetch_pubmed_abstracts(query, max_results=max_results)

    all_points = []

    for paper in papers:

        if not paper["abstract"]:
            continue

        chunks = simple_chunk(paper["abstract"])
        vectors = embed_texts(chunks)

        for i, (chunk, vector) in enumerate(zip(chunks, vectors)):

            payload = build_payload(
                text=chunk,
                pmid=paper["pmid"],
                title=paper["title"],
                journal=paper["journal"],
                year=paper["year"],
                authors=paper["authors"],
                section="Abstract",
                chunk_index=i,
                api_query=query
            )

            point = PointStruct(
                id=str(uuid.uuid4()),
                vector=vector,
                payload=payload
            )

            all_points.append(point)

    if all_points:
        client.upsert(
            collection_name=COLLECTION,
            points=all_points
        )

    print(f"Ingested {len(all_points)} chunks from PubMed.")


if __name__ == "__main__":

    ingest_from_pubmed(
        query="metformin chronic kidney disease",
        max_results=5
    )
