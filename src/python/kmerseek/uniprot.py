"""
Parse JSON file from UniProt
"""


def get_domains(uniprot_id, start, end):
    """For the query uniprot ID (primary accession, e.g. O43236) provided, find overlapping regions

    Probably would do something like this:
    https://github.com/olgabot/botryllus-mhc/blob/main/notebooks/400_query_uniprot_rest_api.ipynb
    Except not querying the REST API because get throttled with too many requests
    Better to pre-download the JSON of the sequence annotation data:
    https://www.uniprot.org/help/sequence_annotation

    Or maybe just download the UniProt/SwissProt XML every 8 weeks from here?
    https://www.uniprot.org/help/downloads
    """

    pass
