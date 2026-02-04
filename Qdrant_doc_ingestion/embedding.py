from sentence_transformers import SentenceTransformer

_model = SentenceTransformer("all-mpnet-base-v2")

def embed_texts(texts):
    return _model.encode(texts).tolist()
