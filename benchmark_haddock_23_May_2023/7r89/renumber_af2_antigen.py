import os

def get_resids(a_file):
    resids = []
    prev_resid = "blu"
    with open(a_file, "r") as rfile:
        for ln in rfile:
            if ln.startswith("ATOM"):
                splt_ln = ln.split()
                resid = f"{splt_ln[4]}-{splt_ln[5]}"
                if resid !=prev_resid:
                    resids.append(resid)
                    prev_resid = resid
    return resids

ant_resids = get_resids("7r89_antigen.pdb")
ant_af2_resids = get_resids("7r89_AF2_antigen_model_orig.pdb")

k = -1
prev_resid = "blu"

with open("7r89_AF2_antigen_model.pdb", "w") as wfile:
    with open("7r89_AF2_antigen_model_orig.pdb", "r") as rfile:
        for ln in rfile:
            if ln.startswith("ATOM"):
                splt_ln = ln.split()
                resid = f"{splt_ln[4]}-{splt_ln[5]}"
                print(resid)
                if resid != prev_resid:
                    prev_resid = resid
                    k+=1
                new_chain = ant_resids[k].split("-")[0]
                new_resid =  ant_resids[k].split("-")[1]
                space_res = 4 - len(new_resid)
                new_ln = f"{ln[:21]}{new_chain}{' ' * space_res}{new_resid}{ln[26:]}"
                wfile.write(new_ln)
            else:
                wfile.write(ln)
