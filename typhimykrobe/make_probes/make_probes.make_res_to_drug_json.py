import json

data = {}

with open("qrdr.acrB.aa.in.txt") as f:
    for line in f:
        gene, var, _ = line.split()
        assert var not in data
        if "acrB" in gene:
            data[var] = ["azithromycin"]
        else:
            data[var] = ["ciprofloxacin"]

#with open("amr_genes_drugs_dfrTest.txt") as f:
with open("amr_genes_drugs.txt") as f:
    for line in f:
        fields = line.strip().split()
        gene = fields[0]
        drug = fields[1]
        if not gene in data.keys():
            data[gene] = [drug]
        else:
            data[gene].append(drug)

print(json.dumps(data, sort_keys=True, indent=2))
