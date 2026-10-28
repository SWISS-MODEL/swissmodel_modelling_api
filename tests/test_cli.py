"""Tests for the command line of the wrapper, up to talking to the API."""

import sys
from unittest import mock

import pytest

import sm_api_wrapper as w

TOKEN_MSG = "SWISS-MODEL token is not provided in credentials!"
CRN_TRG = "KSCCPTTAARNQYNICRLPGTPRPVCAALSGCKIISGTGCPPGYRH"
CRN_TPL = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
ALIGNMENT_OPTS = ["-p", "1crn", "-a", "1", "-c", "A", "-o", "0"]


def run(argv, monkeypatch, capsys):
    """Run the wrapper without token, return exit code and stderr."""
    monkeypatch.delenv("SWISSMODEL_API_TOKEN", raising=False)
    monkeypatch.setattr(sys, "argv", ["sm_api_wrapper.py", *argv])
    with pytest.raises(SystemExit) as exc:
        w._main()
    return exc.value.code, capsys.readouterr().err


@pytest.fixture(name="files")
def fixture_files(tmp_path):
    """Input files for the wrapper, keyed by short names."""
    contents = {
        "aln": f">THN_DENCL\n{CRN_TRG}\n>1crn.1.A offset=0\n{CRN_TPL}\n",
        "aln3": f">t\n{CRN_TRG}\n>1crn.1.A\n{CRN_TPL}\n>extra\nGG\n",
        "fas": ">a\nMKV\n>b\nGGA\n",
        "headers_only": ">a\n>b\n",
        "missing_chain": ">A\nMVVKAVCVINGDAKG\n>B\n",
        "empty": "",
        "pdb": "ATOM      1  N   THR A   1      17.047  14.099   3.625\n",
    }
    paths = {"out": str(tmp_path)}
    for name, text in contents.items():
        path = tmp_path / f"{name}.txt"
        path.write_text(text)
        paths[name] = str(path)
    return paths


@pytest.mark.parametrize(
    "argv",
    [
        ["single", "automodel", "{out}", "MKV", "__gt__b GGA"],
        [
            "single",
            "alignment",
            "{out}",
            "--target-sequence=KSCC",
            "--template-sequence=TTCC",
            *ALIGNMENT_OPTS,
        ],
        [
            "single",
            "alignment",
            "{out}",
            "--target-sequence=-KSCC",
            "--template-sequence=-TTCC",
            *ALIGNMENT_OPTS,
        ],
        ["single", "usertemplate", "{out}", "MKV", "GGA", "-f", "{pdb}"],
        ["batch", "automodel", "{out}", "{fas}"],
        ["batch", "alignment", "{out}", "{aln}"],
        ["batch", "usertemplate", "{out}", "{fas}", "-f", "{pdb}"],
    ],
)
def test_valid_input_reaches_token_check(argv, files, monkeypatch, capsys):
    code, err = run([a.format(**files) for a in argv], monkeypatch, capsys)
    assert code == 1
    assert TOKEN_MSG in err


@pytest.mark.parametrize(
    "argv, message",
    [
        (["single"], "the following arguments are required"),
        (
            ["batch", "automodel", "{out}", "{fas}", "-p", "1crn"],
            "unrecognized arguments: -p 1crn",
        ),
        (
            ["single", "usertemplate", "{out}", "MKV"],
            "the following arguments are required: -f/--template-file",
        ),
        (
            [
                "single",
                "alignment",
                "{out}",
                "--target-sequence=KS",
                *ALIGNMENT_OPTS,
            ],
            "the following arguments are required: -t/--template-sequence",
        ),
        (
            ["single", "automodel", "{out}", "MKV", "--project-title="],
            "Option '--project-title' can not have an empty string as value",
        ),
        (
            ["single", "automodel", "{out}", "__gt__seq1 "],
            "No target sequence given",
        ),
        (
            [
                "single",
                "alignment",
                "{out}",
                "--target-sequence=KS",
                "--template-sequence=__gt__x ",
                *ALIGNMENT_OPTS,
            ],
            "No template sequence given",
        ),
        (
            ["batch", "automodel", "{out}", "{headers_only}"],
            "Target FASTA record 'a' has no sequence",
        ),
        (
            ["batch", "automodel", "{out}", "{missing_chain}"],
            "Target FASTA record 'B' has no sequence",
        ),
        (
            [
                "batch",
                "usertemplate",
                "{out}",
                "{missing_chain}",
                "-f",
                "{pdb}",
            ],
            "Target FASTA record 'B' has no sequence",
        ),
        (
            ["batch", "automodel", "{out}", "{empty}"],
            "No target sequence given",
        ),
        (
            ["batch", "alignment", "{out}", "{aln3}"],
            "Alignment file must contain exactly 2 sequences, found 3",
        ),
    ],
)
def test_invalid_input_exits_2_before_token(
    argv, message, files, monkeypatch, capsys
):
    code, err = run([a.format(**files) for a in argv], monkeypatch, capsys)
    assert code == 2
    assert message in err
    assert TOKEN_MSG not in err


def test_batch_alignment_payload(tmp_path, monkeypatch, capsys):
    aln = tmp_path / "aln.fas"
    aln.write_text(f">t\n{CRN_TRG}\n>1CRN.1.A offset=3\n{CRN_TPL}\n")
    submitted = mock.Mock(ok=True, status_code=202)
    submitted.json.return_value = {"project_id": "p1"}
    details = mock.Mock(ok=True, status_code=200)
    details.json.return_value = {
        "project_id": "p1",
        "target_sequences": [CRN_TRG],
        "status": "COMPLETED",
        "message": "",
        "date_created": "2026-10-28T11:00:00",
        "models": [],
        "project_title": None,
        "view_url": "https://swissmodel.expasy.org/interactive/p1/",
    }
    post = mock.Mock(return_value=submitted)
    monkeypatch.setattr(w.requests, "post", post)
    monkeypatch.setattr(w.requests, "get", mock.Mock(return_value=details))
    monkeypatch.setattr(w.time, "sleep", lambda _: None)
    monkeypatch.setenv("SWISSMODEL_API_TOKEN", "token")
    monkeypatch.setattr(
        sys,
        "argv",
        ["sm_api_wrapper.py", "batch", "alignment", str(tmp_path), str(aln)],
    )

    with pytest.raises(SystemExit) as exc:
        w._main()

    assert exc.value.code == 0
    assert post.call_args.kwargs["json"] == {
        "assembly_id": 1,
        "auth_asym_id": "A",
        "pdb_id": "1crn",
        "target_sequences": [CRN_TRG],
        "template_seqres_offset": 3,
        "template_sequence": CRN_TPL,
        "project_title": None,
    }
