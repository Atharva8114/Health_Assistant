import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import time
from Bio import Entrez
from bs4 import BeautifulSoup

# ALWAYS use your email so NCBI doesn't block you
Entrez.email = "your_email@example.com" 

def search_pmc_articles(query, max_results=5):
    """
    Searches PubMed Central (PMC) for OPEN ACCESS full-text articles.
    """
    print(f"🔎 Searching PMC for: '{query}'...")
    
    # We add "open access[filter]" to ensure we can actually read the full text
    full_query = f"{query} AND open access[filter]"
    
    try:
        # 1. Search for IDs
        handle = Entrez.esearch(db="pmc", term=full_query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        print(f"✅ Found {len(id_list)} full-text articles.")
        return id_list
        
    except Exception as e:
        print(f"❌ Error searching PMC: {e}")
        return []

def fetch_pmc_details(pmc_id):
    """
    Downloads the Full Text XML and cleans it.
    """
    try:
        # 2. Fetch the Full Record (XML)
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="xml")
        xml_data = handle.read()
        handle.close()
        
        # 3. Parse XML with BeautifulSoup
        soup = BeautifulSoup(xml_data, "xml")
        
        # --- Metadata Extraction ---
        article = {}
        article['pmid'] = pmc_id # Using PMC ID as the identifier
        
        # Get Title
        title_tag = soup.find("article-title")
        article['title'] = title_tag.get_text() if title_tag else "No Title"
        
        # Get Journal
        journal_tag = soup.find("journal-title")
        article['journal'] = journal_tag.get_text() if journal_tag else "Unknown Journal"
        
        # Get Year
        pub_date = soup.find("pub-date", {"pub-type": "epub"})
        if not pub_date:
            pub_date = soup.find("pub-date")
        
        if pub_date and pub_date.find("year"):
            article['year'] = int(pub_date.find("year").get_text())
        else:
            article['year'] = 0

        # --- CRITICAL: FULL TEXT EXTRACTION ---
        # PMC stores the main text in the <body> tag.
        body = soup.find("body")
        
        if body:
            # We iterate through sections to make it readable
            # Removing citations like [1] or (Smith et al) would happen here in a pro app,
            # but for now, we just get the raw text paragraphs.
            text_sections = []
            
            for sec in body.find_all("sec"):
                title = sec.find("title")
                section_title = title.get_text() if title else "Section"
                
                # Get all paragraphs in this section
                paragraphs = [p.get_text() for p in sec.find_all("p")]
                if paragraphs:
                    # Format: "## Introduction \n Text..."
                    text_sections.append(f"## {section_title}\n" + "\n".join(paragraphs))
            
            full_text = "\n\n".join(text_sections)
            article['abstract'] = full_text # We save full text in the 'abstract' field to keep compatibility
            
        else:
            # Fallback if no body found (some older docs)
            article['abstract'] = "Full text not available in XML format."

        article['study_type'] = "Full Text Article"
        
        return article

    except Exception as e:
        print(f"⚠️ Error processing PMC{pmc_id}: {e}")
        return None

def fetch_all_pmc_articles(query="metformin", max_results=5):
    pmc_ids = search_pmc_articles(query, max_results)
    articles = []
    
    for pmc_id in pmc_ids:
        data = fetch_pmc_details(pmc_id)
        if data and len(data['abstract']) > 500: # Only keep if we actually got text
            articles.append(data)
            print(f"   📄 Processed: {data['title'][:50]}...")
        time.sleep(0.5) # Be polite to the API
        
    return articles

# --- TEST BLOCK ---
if __name__ == "__main__":
    results = fetch_all_pmc_articles("metformin side effects", 1)
    if results:
        print("\n--- SAMPLE EXTRACTED TEXT ---")
        print(results[0]['abstract'][:1000]) # Print first 1000 characters

