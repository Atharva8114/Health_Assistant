from gliner import GLiNER
import torch

# Load a lightweight, high-performance model
# "urchade/gliner_small-v2.1" is fast and accurate for this task
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"[INFO] Loading NER model on {device}...")
model = GLiNER.from_pretrained("urchade/gliner_small-v2.1").to(device)

def extract_medical_entities(text: str):
    """
    Scans text and returns a dictionary of lists: 
    {'drugs': [], 'conditions': [], 'biomarkers': [], 'symptoms': []}
    """
    if not text:
        return {"drugs": [], "conditions": [], "biomarkers": [], "symptoms": []}

    # Define what we are looking for
    labels = ["drug", "medical condition", "biomarker", "symptom"]
    
    # Predict
    entities = model.predict_entities(text, labels, threshold=0.3)
    
    # Organize results
    results = {
        "drugs": [],
        "conditions": [],
        "biomarkers": [],
        "symptoms": []
    }
    
    for entity in entities:
        text_clean = entity["text"].lower().strip()
        label = entity["label"]
        
        if label == "drug" and text_clean not in results["drugs"]:
            results["drugs"].append(text_clean)
        elif label == "medical condition" and text_clean not in results["conditions"]:
            results["conditions"].append(text_clean)
        elif label == "biomarker" and text_clean not in results["biomarkers"]:
            results["biomarkers"].append(text_clean)
        elif label == "symptom" and text_clean not in results["symptoms"]:
            results["symptoms"].append(text_clean)
            
    return results