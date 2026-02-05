from langchain_text_splitters import RecursiveCharacterTextSplitter

# CHANGED: chunk_size=1000 (Fits inside NER model limit)
def simple_chunk(text, chunk_size=1000, overlap=100):
    splitter = RecursiveCharacterTextSplitter(
        chunk_size=chunk_size,
        chunk_overlap=overlap,
        separators=["\n\n", "\n", ".", " ", ""]
    )
    return splitter.split_text(text)