from click.testing import CliRunner

from kmerseek.main import cli


def test_search(ced9, bcl2_first25, bcl2_hp_k16_sig_zip):
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(
        cli,
        [
            "search",
            # "--output",
            # "kmerseek_results.csv",
            "--ksize",
            "16",
            ced9,
            bcl2_first25,
        ],
    )
    assert result.exit_code == 0

    assert (
        result.stdout
        == """match_name,query_name,query_start,query_end,query,match_start,match_end,match
sp|Q9UK96|FBX10_HUMAN F-box only protein 10 OS=Homo sapiens OX=9606 GN=FBXO10 PE=1 SV=3,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,59,92,MSIGESIDGKINDWEEPGIVGVVVCGRMMFSLK,57,90,PNWPNQPDVEPESWREAAGIYILYHGNPVVSGN
sp|Q12982|BNIP2_HUMAN BCL2/adenovirus E1B 19 kDa protein-interacting protein 2 OS=Homo sapiens OX=9606 GN=BNIP2 PE=1 SV=1,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,76,108,RLDIEGFVVDYFTHRILFVYTSLFIKTRIRNN,23,55,SIEADILAITGPEDQPLLAVTRPFISSKFSQK
"sp|Q9BXH1|BBC3_HUMAN Bcl-2-binding component 3, isoforms 1/2 OS=Homo sapiens OX=9606 GN=BBC3 PE=1 SV=1",sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,170,187,LIGLISFGGFVAAKMME,46,63,APAAPTLLPAAYLCAPT
sp|Q13625|ASPP2_HUMAN Apoptosis-stimulating of p53 protein 2 OS=Homo sapiens OX=9606 GN=TP53BP2 PE=1 SV=2,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,241,257,KVGRRKQNRRWSMIGA,1084,1100,TIIHREDEDEIEWWWA
sp|Q16611|BAK_HUMAN Bcl-2 homologous antagonist/killer OS=Homo sapiens OX=9606 GN=BAK1 PE=1 SV=1,sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1,245,261,RKQNRRWSMIGAGVTA,42,58,HQQEQEAEGVAAPADP
"""
    )
    assert (
        result.stderr
        == """
---
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
match: HQQEQEAEGVAAPADP (42-58)

"""
    )
