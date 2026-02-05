import uuid
from typing import List, Dict, Any

from qdrant_client.models import PointStruct
from qdrant_client.http.exceptions import UnexpectedResponse

from pubmed_fetcher import fetch_all_pmc_articles
from chunker import simple_chunk
from embedding import embed_texts
from schema_builder import build_payload
# NEW IMPORT (Ensure entity_extractor.py exists)
from entity_extractor import extract_medical_entities 
from qdrant_store import (
    get_client,
    create_collection_if_not_exists,
    create_indexes,
    COLLECTION
)

def ingest_from_pubmed(query: str, max_results: int = 5) -> None:
    """
    Orchestrates the fetching, processing, and ingestion of PMC full-text papers.
    """
    try:
        client = get_client()
        create_collection_if_not_exists(client)
        create_indexes(client)
    except Exception as e:
        print(f"[ERROR] Failed to initialize Qdrant connection: {e}")
        return

    print(f"[INFO] Querying PMC for full-text articles: '{query}'")
    
    try:
        papers = fetch_all_pmc_articles(query, max_results=max_results)
    except Exception as e:
        print(f"[ERROR] Failed to fetch articles from PMC: {e}")
        return

    if not papers:
        print("[WARNING] No papers found matching the query.")
        return

    all_points = []
    
    for paper in papers:
        pmid = paper.get("pmid", "unknown")
        title = paper.get("title", "No Title")
        abstract_text = paper.get("abstract")

        if not abstract_text:
            print(f"[WARNING] Skipping paper {pmid}: No text content available.")
            continue

        print(f"[INFO] Processing paper: {title[:50]}...")

        # Chunking
        chunks = simple_chunk(abstract_text)
        if not chunks:
            print(f"[WARNING] Skipping paper {pmid}: Text chunking returned empty list.")
            continue

        # Embedding
        try:
            vectors = embed_texts(chunks)
        except Exception as e:
            print(f"[ERROR] Embedding failed for paper {pmid}: {e}")
            continue

        if len(chunks) != len(vectors):
            print(f"[ERROR] Mismatch in chunks and vectors count for paper {pmid}. Skipping.")
            continue

        # Point Construction
        for i, (chunk, vector) in enumerate(zip(chunks, vectors)):
            try:
                # --- NEW STEP: Extract Medical Entities ---
                # We extract specific terms (Drugs, Symptoms) from this specific chunk
                extracted_data = extract_medical_entities(chunk)
                # ------------------------------------------

                payload = build_payload(
                    text=chunk,
                    pmid=pmid,
                    title=title,
                    journal=paper.get("journal", "Unknown"),
                    year=paper.get("year", 0),
                    authors=paper.get("authors", ["PMC Full Text"]),
                    section="Full Text",
                    chunk_index=i,
                    api_query=query,
                    entities=extracted_data  # <--- PASSING DATA HERE
                )

                point = PointStruct(
                    id=str(uuid.uuid4()),
                    vector=vector,
                    payload=payload
                )
                all_points.append(point)
            except Exception as e:
                print(f"[ERROR] Failed to build payload for chunk {i} of paper {pmid}: {e}")

    # Batch Upload
    if all_points:
        try:
            print(f"[INFO] Uploading {len(all_points)} vectors to collection '{COLLECTION}'...")
            client.upsert(
                collection_name=COLLECTION,
                points=all_points
            )
            print(f"[SUCCESS] Ingestion complete. Total chunks indexed: {len(all_points)}")
        except UnexpectedResponse as e:
            print(f"[ERROR] Qdrant upsert failed: {e}")
        except Exception as e:
            print(f"[ERROR] An unexpected error occurred during upload: {e}")
    else:
        print("[WARNING] No valid data points processed to upload.")

if __name__ == "__main__":
    ingest_from_pubmed(
        query="metformin side effects",
        max_results=2 # Kept low for testing speed
    )

