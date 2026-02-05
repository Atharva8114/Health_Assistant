from sentence_transformers import SentenceTransformer

_model = SentenceTransformer("BAAI/bge-m3")

def embed_texts(texts):
    return _model.encode(texts, normalize_embeddings=True).tolist()
