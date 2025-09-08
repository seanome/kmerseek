import os
import polars as pl
from io import StringIO
from polars.testing import assert_frame_equal
import subprocess
import sys


def test_search(ced9, bcl2_first25, bcl2_hp_k16_sig_zip):
    # Use subprocess to properly capture stderr separately
    cmd = [
        sys.executable,
        "-m",
        "kmerseek.main",
        "search",
        "--ksize",
        "16",
        ced9,
        bcl2_first25,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd())
    assert result.returncode == 0

    # TODO: Actually test for the signature in this file by reading it in with
    # `sourmash.load_file_as_signatures`
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{ced9}.hp.k16.scaled5.sig.zip")
    assert os.path.exists(f"{bcl2_first25}.hp.k16.scaled5.sig.zip")

    true_output = pl.read_csv(
        StringIO(
            """query_name,query_md5,match_name,containment,intersect_hashes,ksize,scaled,moltype,match_md5,jaccard,max_containment,average_abund,median_abund,std_abund,query_containment_ani,match_containment_ani,average_containment_ani,max_containment_ani,n_weighted_found,total_weighted_hashes,containment_target_in_query,f_weighted_target_in_query
sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,fe3714626e8180caf90f78091563aae6,sp|Q12982|BNIP2_HUMAN BCL2/adenovirus E1B 19 kDa protein-interacting protein 2 OS=Homo sapiens OX=9606 GN=BNIP2 PE=1 SV=1,0.04081632653061224,2,48,5,hp,7bbc6e2ea3a472034fc31321943032ee,0.02040816326530612,0.04081632653061224,1.0,1.0,0.0,0.9355328459682174,0.934753456124389,0.9351431510463032,0.9355328459682174,2,51,0.0392156862745098,0.0392156862745098
sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,fe3714626e8180caf90f78091563aae6,sp|Q13625|ASPP2_HUMAN Apoptosis-stimulating of p53 protein 2 OS=Homo sapiens OX=9606 GN=TP53BP2 PE=1 SV=2,0.02040816326530612,1,48,5,hp,35da5dcf3561c6c0b0aaa34a118eabef,0.0036101083032490976,0.02040816326530612,1.0,1.0,0.0,0.9221202973899911,0.8929697781452893,0.9075450377676402,0.9221202973899911,1,230,0.004366812227074236,0.004347826086956522
sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,fe3714626e8180caf90f78091563aae6,sp|Q16611|BAK_HUMAN Bcl-2 homologous antagonist/killer OS=Homo sapiens OX=9606 GN=BAK1 PE=1 SV=1,0.02040816326530612,1,48,5,hp,1f59cdb10b02a7c6baff18b034518599,0.011111111111111112,0.023809523809523808,1.0,1.0,0.0,0.9221202973899911,0.9250864216273635,0.9236033595086773,0.9250864216273635,1,42,0.023809523809523808,0.023809523809523808
sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,fe3714626e8180caf90f78091563aae6,"sp|Q9BXH1|BBC3_HUMAN Bcl-2-binding component 3, isoforms 1/2 OS=Homo sapiens OX=9606 GN=BBC3 PE=1 SV=1",0.04081632653061224,2,48,5,hp,1d49aa1205276b9ba0176c6680cacd6d,0.024390243902439025,0.05714285714285714,1.0,1.0,0.0,0.9355328459682174,0.9421138187376149,0.9388233323529162,0.9421138187376149,2,35,0.05714285714285714,0.05714285714285714
sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,fe3714626e8180caf90f78091563aae6,sp|Q9UK96|FBX10_HUMAN F-box only protein 10 OS=Homo sapiens OX=9606 GN=FBXO10 PE=1 SV=3,0.061224489795918366,3,48,5,hp,97f5f83c6214d6792113785b96747383,0.014354066985645933,0.061224489795918366,1.0,1.0,0.0,0.9434689410983454,0.9201376138657374,0.9318032774820415,0.9434689410983454,3,164,0.018404907975460124,0.018292682926829267
"""
        )
    ).sort("match_name")

    # Extract CSV data from stdout (skip any non-CSV content that comes before it)
    stdout_lines = result.stdout.split("\n")
    csv_start = None
    for i, line in enumerate(stdout_lines):
        if line.startswith("query_name,query_md5,match_name"):
            csv_start = i
            break

    if csv_start is not None:
        csv_content = "\n".join(stdout_lines[csv_start:])
        test_output = pl.read_csv(StringIO(csv_content)).sort("match_name")
    else:
        # Fallback: try to parse the entire stdout with truncate_ragged_lines
        test_output = pl.read_csv(
            StringIO(result.stdout), truncate_ragged_lines=True
        ).sort("match_name")

    assert_frame_equal(test_output, true_output)


def test_search_extract_kmers(ced9, bcl2_first25, bcl2_hp_k16_sig_zip):
    # Use subprocess to properly capture stderr separately
    cmd = [
        sys.executable,
        "-m",
        "kmerseek.main",
        "search",
        "--extract-kmers",
        "--ksize",
        "16",
        ced9,
        bcl2_first25,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd())
    assert result.returncode == 0

    # TODO: Actually test for the signature in this file by reading it in with
    # `sourmash.load_file_as_signatures`
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{ced9}.hp.k16.scaled5.sig.zip")
    assert os.path.exists(f"{bcl2_first25}.hp.k16.scaled5.sig.zip")

    true_output = pl.read_csv(
        StringIO(
            """match_name,query_name,query_start,query_end,query,match_start,match_end,match,encoded,length
sp|Q12982|BNIP2_HUMAN BCL2/adenovirus E1B 19 kDa protein-interacting protein 2 OS=Homo sapiens OX=9606 GN=BNIP2 PE=1 SV=1,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,76,108,RLDIEGFVVDYFTHRILFVYTSLFIKTRIRNN,23,55,SIEADILAITGPEDQPLLAVTRPFISSKFSQK,phphphhhhphhppphhhhhpphhhppphppp,32
sp|Q13625|ASPP2_HUMAN Apoptosis-stimulating of p53 protein 2 OS=Homo sapiens OX=9606 GN=TP53BP2 PE=1 SV=2,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,241,257,KVGRRKQNRRWSMIGA,1084,1100,TIIHREDEDEIEWWWA,phhppppppphphhhh,16
sp|Q16611|BAK_HUMAN Bcl-2 homologous antagonist/killer OS=Homo sapiens OX=9606 GN=BAK1 PE=1 SV=1,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,245,261,RKQNRRWSMIGAGVTA,42,58,HQQEQEAEGVAAPADP,pppppphphhhhhhph,16
"sp|Q9BXH1|BBC3_HUMAN Bcl-2-binding component 3, isoforms 1/2 OS=Homo sapiens OX=9606 GN=BBC3 PE=1 SV=1",sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,170,187,LIGLISFGGFVAAKMME,46,63,APAAPTLLPAAYLCAPT,hhhhhphhhhhhhphhp,17
sp|Q9UK96|FBX10_HUMAN F-box only protein 10 OS=Homo sapiens OX=9606 GN=FBXO10 PE=1 SV=3,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,59,92,MSIGESIDGKINDWEEPGIVGVVVCGRMMFSLK,57,90,PNWPNQPDVEPESWREAAGIYILYHGNPVVSGN,hphhpphphphpphpphhhhhhhhphphhhphp,33
"""
        )
    ).sort("match_name")

    # CSV data should be in stdout, visual output should be in stderr
    test_output = pl.read_csv(StringIO(result.stdout)).sort("match_name")

    assert_frame_equal(test_output, true_output)

    assert (
        """---
Query Name: sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
Match Name: sp|Q9UK96|FBX10_HUMAN F-box only protein 10 OS=Homo sapiens OX=9606 GN=FBXO10 PE=1 SV=3
query: MSIGESIDGKINDWEEPGIVGVVVCGRMMFSLK (59-92)
alpha: hphhpphphphpphpphhhhhhhhphphhhphp
match: PNWPNQPDVEPESWREAAGIYILYHGNPVVSGN (57-90)

---
Query Name: sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
Match Name: sp|Q12982|BNIP2_HUMAN BCL2/adenovirus E1B 19 kDa protein-interacting protein 2 OS=Homo sapiens OX=9606 GN=BNIP2 PE=1 SV=1
query: RLDIEGFVVDYFTHRILFVYTSLFIKTRIRNN (76-108)
alpha: phphphhhhphhppphhhhhpphhhppphppp
match: SIEADILAITGPEDQPLLAVTRPFISSKFSQK (23-55)

---
Query Name: sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
Match Name: sp|Q9BXH1|BBC3_HUMAN Bcl-2-binding component 3, isoforms 1/2 OS=Homo sapiens OX=9606 GN=BBC3 PE=1 SV=1
query: LIGLISFGGFVAAKMME (170-187)
alpha: hhhhhphhhhhhhphhp
match: APAAPTLLPAAYLCAPT (46-63)

---
Query Name: sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
Match Name: sp|Q13625|ASPP2_HUMAN Apoptosis-stimulating of p53 protein 2 OS=Homo sapiens OX=9606 GN=TP53BP2 PE=1 SV=2
query: KVGRRKQNRRWSMIGA (241-257)
alpha: phhppppppphphhhh
match: TIIHREDEDEIEWWWA (1084-1100)

---
Query Name: sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
Match Name: sp|Q16611|BAK_HUMAN Bcl-2 homologous antagonist/killer OS=Homo sapiens OX=9606 GN=BAK1 PE=1 SV=1
query: RKQNRRWSMIGAGVTA (245-261)
alpha: pppppphphhhhhhph
match: HQQEQEAEGVAAPADP (42-58)"""
        in result.stderr
    )
