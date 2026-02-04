from qdrant_client import QdrantClient
from qdrant_client.models import VectorParams, Distance, PayloadSchemaType

COLLECTION = "medical_documents"


def get_client():
    return QdrantClient(host="localhost", port=6333)


def create_collection_if_not_exists(client):
    collections = client.get_collections().collections
    names = [c.name for c in collections]

    if COLLECTION not in names:
        client.create_collection(
            collection_name=COLLECTION,
            vectors_config=VectorParams(
                size=768,
                distance=Distance.COSINE
            )
        )


def create_indexes(client):
    client.create_payload_index(
        collection_name=COLLECTION,
        field_name="pmid",
        field_schema=PayloadSchemaType.KEYWORD
    )

    client.create_payload_index(
        collection_name=COLLECTION,
        field_name="year",
        field_schema=PayloadSchemaType.INTEGER
    )

    client.create_payload_index(
        collection_name=COLLECTION,
        field_name="journal",
        field_schema=PayloadSchemaType.KEYWORD
    )

    client.create_payload_index(
        collection_name=COLLECTION,
        field_name="study_type",
        field_schema=PayloadSchemaType.KEYWORD
    )

    client.create_payload_index(
        collection_name=COLLECTION,
        field_name="entities.drugs",
        field_schema=PayloadSchemaType.KEYWORD
    )
