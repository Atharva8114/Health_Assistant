import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from Bio import Entrez
import xml.etree.ElementTree as ET


# NCBI requires an email
Entrez.email = "just_put_some_temp_mail"


def fetch_pubmed_abstracts(query, max_results=5):
    """
    Returns list of dicts:
    [
      {
        pmid, title, abstract, journal, year, authors
      }
    ]
    """

    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=max_results
    )

    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]

    if not id_list:
        return []

    fetch_handle = Entrez.efetch(
        db="pubmed",
        id=",".join(id_list),
        rettype="xml",
        retmode="text"
    )

    data = fetch_handle.read()
    fetch_handle.close()

    root = ET.fromstring(data)

    papers = []

    for article in root.findall(".//PubmedArticle"):

        pmid = article.findtext(".//PMID")

        title = article.findtext(".//ArticleTitle")

        abstract_parts = article.findall(".//AbstractText")
        abstract = " ".join([a.text for a in abstract_parts if a.text])

        journal = article.findtext(".//Journal/Title")

        year = article.findtext(".//PubDate/Year")

        authors = []
        for a in article.findall(".//Author"):
            last = a.findtext("LastName")
            fore = a.findtext("ForeName")
            if last and fore:
                authors.append(f"{fore} {last}")

        papers.append({
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": int(year) if year and year.isdigit() else None,
            "authors": authors
        })

    return papers
