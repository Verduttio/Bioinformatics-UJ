from Bio import Entrez, SeqIO


#### ZADANIE 1 ####
def get_gene_ids_connected_with_brca1_homo_sapiens():
    Entrez.email = "bszwaja20@gmail.com"

    query = "BRCA1[Gene] AND Homo sapiens[Organism]"
    handle = Entrez.esearch(db="gene", term=query)
    record = Entrez.read(handle)

    return record["IdList"]


def get_gene_record(gene_id):
    handle = Entrez.esummary(db="gene", id=gene_id, rettype="gb", retmode="text")
    gene_record = Entrez.read(handle)
    handle.close()

    return gene_record["DocumentSummarySet"]["DocumentSummary"][0]


def show_gene_record_details(gene_record):
    print("Name: ", gene_record["Name"])
    print("Chromosome: ", gene_record["Chromosome"])
    print("MapLocation: ", gene_record["MapLocation"])
    print("Description: ", gene_record["Description"])


#### ZADANIE 3 ####
def get_omim_ids_connected_with_gene_id(gene_id):
    handle = Entrez.elink(dbfrom="gene", db="omim", id=gene_id)
    result = Entrez.read(handle)
    omim_ids = []
    for link in result[0]["LinkSetDb"][0]["Link"]:
        omim_ids.append(link["Id"])

    return omim_ids


def get_omim_records_from_omim_ids(omim_ids):
    handle = Entrez.esummary(db="omim", id=",".join(omim_ids))
    return Entrez.read(handle)


#### ZADANIE 4 ####
def get_protein_ids_connected_with_gene_id(gene_id):
    handle = Entrez.elink(dbfrom="gene", db="protein", id=gene_id, retmax=10)
    result = Entrez.read(handle)
    protein_ids = []
    for link in result[0]["LinkSetDb"][0]["Link"]:
        protein_ids.append(link["Id"])

    return protein_ids


def get_gi_ids_from_protein_ids(protein_ids):
    handle = Entrez.esummary(db="protein", id=",".join(protein_ids), retmax=10)
    gi_ids = []
    for record in Entrez.read(handle):
        gi_ids.append(int(record["Gi"]))

    return gi_ids


#### ZADANIE 5 ####
def get_protein_121949022():
    handle = Entrez.efetch(db="protein", id="121949022", rettype="gb", retmode="text")
    return SeqIO.read(handle, "genbank")


#### ZADANIE 6 ####
def get_mutations_of_brca1_for_human(gene_id="672"):
    handle = Entrez.elink(dbfrom="gene", db="snp", id=gene_id, term="Homo sapiens")
    result = Entrez.read(handle)
    mutation_ids = []
    for link in result[0]["LinkSetDb"][0]["Link"]:
        mutation_ids.append(link["Id"])

    return mutation_ids


def print_mutation_info(mutation_ids):
    for mutation_id in mutation_ids:
        handle = Entrez.esummary(db="snp", id=mutation_id)
        record = Entrez.read(handle)
        record = record["DocumentSummarySet"]["DocumentSummary"][0]
        print("SNP_ID: ", record["SNP_ID"])
        print("SNP_CLASS: ", record["SNP_CLASS"])
        print("GENES: ", record["GENES"])
        print("CONTIGPOS: ", record["CHRPOS"])
        print()


#### ZADANIE 7 ####
def conclusions():
    print("""
    Na początku znaleźliśmy gen u człowieka, który jest powiązany z rakiem piersi.
    Następnie na podstawie znalezionego genu wyszukaliśmy inne choroby za które może on również odpowiadać.
    Idąc dalej znaleźliśmy białka związane z tym genem. Było ich dosyć sporo, więc wypisaliśmy tylko 10 pierwszych.
    Następnie znaleźliśmy konkretny białko i wypisaliśmy jego nazwę i sekwencję.
    Na koniec znaleźliśmy mutacje tego genu u człowieka i dla pierwszych 50 wypisaliśmy informacje o nich.
    Było ich bardzo wiele, może to sugerować, że gen ten nie jest bardzo rzaedki i może być odpowiedzialny za niektóre choroby.""")


if __name__ == '__main__':
    gene_ids = get_gene_ids_connected_with_brca1_homo_sapiens()
    gene_records = []
    for gene_id in gene_ids:
        gene_record = get_gene_record(gene_id)
        gene_records.append(gene_record)
        show_gene_record_details(gene_record)

    print()

    # Mamy tylko jeden gen, gdyż wypisały się dane tylko dla jednego genu.
    gene_id = gene_ids[0]
    gene_record = gene_records[0]
    print("Gene id: ", gene_id)

    print()
    #### ZADANIE 3 ####
    omim_ids = get_omim_ids_connected_with_gene_id(gene_id)
    omim_records = get_omim_records_from_omim_ids(omim_ids)

    print("Omim db records for gene id: ", gene_id)
    for omim_record in omim_records:
        print("Title: ", omim_record["Title"])


    print()
    print("First 10 Gi id numbers of protein")
    #### ZADANIE 4 ####
    protein_ids = get_protein_ids_connected_with_gene_id(gene_id)
    gi_ids = get_gi_ids_from_protein_ids(protein_ids)
    print(gi_ids)

    print()
    #### ZADANIE 5 ####
    protein_121949022 = get_protein_121949022()
    print("Protein 121949022")
    print("Name: ", protein_121949022.name)
    print("Sequence: ", protein_121949022.seq)

    print()
    #### ZADANIE 6 ####
    mutation_ids = get_mutations_of_brca1_for_human()
    print("Number of BRCA1 mutations: ", len(mutation_ids))
    print_mutation_info(mutation_ids[:50])

    print()
    #### ZADANIE 7 ####
    conclusions()



