"""Tests for parsing sequences from text fields and FASTA files."""

import re

import pytest

import sm_api_wrapper as w

CRN_TRG = "KSCCPTTAARNQYNICRLPGTPRPVCAALSGCKIISGTGCPPGYRH"
CRN_TPL = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"


def test_defastarise_sequences_2_0_0_behaviour():
    assert w._defastarise_sequences(["MKV"]) == ["MKV"]
    assert w._defastarise_sequences(["__gt__seq1 MKV"]) == ["MKV"]
    assert w._defastarise_sequences([">seq1 MKV LLA"]) == ["MKVLLA"]
    assert w._defastarise_sequences(["MKV LLA"]) == ["MKVLLA"]
    assert w._defastarise_sequences(["MKV\tLL A", "GGA"]) == [
        "MKVLLA",
        "GGA",
    ]
    # known, accepted limitation of the 2.0.0 logic
    assert w._defastarise_sequences([">sp|P1 desc MKV"]) == ["descMKV"]


def test_read_fasta_records(tmp_path):
    fas = tmp_path / "a.fas"
    fas.write_text(">t1 desc\nMK V\nLL\n\n>t2\nGG\n>t3\n")
    assert w._read_fasta(str(fas)) == [
        ("t1 desc", "MKVLL"),
        ("t2", "GG"),
        ("t3", ""),
    ]
    with pytest.raises(
        ValueError, match=re.escape("Target FASTA record 't3' has no sequence")
    ):
        w._defastarise_targets(str(fas))


def test_defastarise_targets(tmp_path):
    fas = tmp_path / "a.fas"
    fas.write_text(">t1 desc\nMK V\nLL\n\n>t2\nGG\n")
    assert w._defastarise_targets(str(fas)) == ["MKVLL", "GG"]


def test_read_fasta_byte_order_mark(tmp_path):
    fas = tmp_path / "bom.fas"
    fas.write_bytes(b"\xef\xbb\xbf>t1\nMKV\n")
    assert w._read_fasta(str(fas)) == [("t1", "MKV")]


def _aln(tmp_path, text, newline="\n"):
    fas = tmp_path / "aln.fas"
    fas.write_bytes(text.replace("\n", newline).encode())
    return str(fas)


def test_parse_alignment_file(tmp_path):
    path = _aln(
        tmp_path,
        f">THN_DENCL\n{CRN_TRG}\n>1crn.1.A offset=12\n{CRN_TPL}\n",
    )
    assert w._parse_alignment_file(path) == {
        "target_sequence": CRN_TRG,
        "template_sequence": CRN_TPL,
        "pdb_id": "1crn",
        "assembly_id": 1,
        "auth_asym_id": "A",
        "template_seqres_offset": 12,
    }


def test_parse_alignment_file_defaults_and_variants(tmp_path):
    # no offset -> 0; gaps kept; upper-case PDB ID kept as given
    path = _aln(tmp_path, ">t\nKS-CC\n>1CRN.2.AB\nTT-CC\n")
    res = w._parse_alignment_file(path)
    assert (
        res["pdb_id"],
        res["assembly_id"],
        res["auth_asym_id"],
        res["template_seqres_offset"],
    ) == ("1CRN", 2, "AB", 0)
    assert res["target_sequence"] == "KS-CC"


def test_parse_alignment_file_description_before_offset(tmp_path):
    path = _aln(tmp_path, ">t\nKSCC\n>1crn.1.A crambin offset=5\nTTCC\n")
    assert w._parse_alignment_file(path)["template_seqres_offset"] == 5


def test_parse_alignment_file_crlf_and_blank_lines(tmp_path):
    path = _aln(tmp_path, ">t\nKSCC\n\n>1crn.1.A\nTTCC\n\n", newline="\r\n")
    res = w._parse_alignment_file(path)
    assert (res["target_sequence"], res["template_sequence"]) == (
        "KSCC",
        "TTCC",
    )


@pytest.mark.parametrize(
    "text, message",
    [
        (
            ">t\nKSCC\n",
            "Alignment file must contain exactly 2 sequences, found 1",
        ),
        (
            ">t\nKS\n>1crn.1.A\nTT\n>x\nGG\n",
            "Alignment file must contain exactly 2 sequences, found 3",
        ),
        (
            ">t\nKSCC\n>1crnA\nTTCC\n",
            "Template header must start with an SMTL ID like '1crn.1.A', "
            + "found '1crnA'",
        ),
        (
            ">t\nKSCC\n>1crn.1.A offset=abc\nTTCC\n",
            "Template offset must be a non-negative integer, found "
            + "'offset=abc'",
        ),
        (">t\n\n>1crn.1.A\nTTCC\n", "Alignment file has an empty sequence"),
        (
            ">t\nKSCC\n>1crn.1.A OFFSET=5\nTTCC\n",
            "Template offset must be written as 'offset=X', found 'OFFSET=5'",
        ),
        (
            ">t\nKSCC\n>1crn.1.A offset = 5\nTTCC\n",
            "Template offset must be written as 'offset=X', found 'offset'",
        ),
    ],
)
def test_parse_alignment_file_invalid(tmp_path, text, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        w._parse_alignment_file(_aln(tmp_path, text))
