from langchain_text_splitters import RecursiveCharacterTextSplitter

def simple_chunk(text, chunk_size=2000, overlap=200):
    splitter = RecursiveCharacterTextSplitter(
        chunk_size=chunk_size,
        chunk_overlap=overlap,
        separators=["\n\n", "\n", ".", " ", ""]
    )
    # Returns a list of strings, just like your old function did
    return splitter.split_text(text)
