from embedding import embed_texts
import numpy as np

texts = [
    "Metformin is used to treat type 2 diabetes.",
    "Metformin helps control blood sugar levels.",
    "Insulin therapy is used for diabetes.",
    "The Eiffel Tower is located in Paris."
]

vectors = embed_texts(texts)

print("✅ Number of embeddings:", len(vectors))
print("📐 Vector dimension:", len(vectors[0]))

# Cosine similarity
def cosine_sim(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

pairs = [
    (0, 1, "Metformin vs Metformin-related"),
    (0, 2, "Metformin vs Insulin"),
    (0, 3, "Metformin vs Eiffel Tower"),
]

print("\n🔍 Similarity checks:")
for i, j, label in pairs:
    sim = cosine_sim(vectors[i], vectors[j])
    print(f"{label}: {sim:.4f}")
