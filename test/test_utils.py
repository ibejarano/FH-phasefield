import json
import pytest
from utils import read_data

def test_read_data_basic(tmp_path, monkeypatch):
    # Create a temporary JSON file in the expected data/ directory
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    fname = "testfile"
    file_path = data_dir / f"{fname}.json"
    content = {"a": 1, "b": 2}
    with open(file_path, "w") as f:
        json.dump(content, f)

    # Patch working directory so read_data finds the file
    monkeypatch.chdir(tmp_path)


    result = read_data(fname)
    assert result == content

def test_read_data_with_override(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    fname = "testfile"
    file_path = data_dir / f"{fname}.json"
    content = {"a": 1, "b": 2}
    with open(file_path, "w") as f:
        json.dump(content, f)

    monkeypatch.chdir(tmp_path)

    result = read_data(fname, overrrides=["a=42"])
    assert result["a"] == 42
    assert result["b"] == 2

def test_read_data_override_nonexistent_key(tmp_path, monkeypatch, capsys):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    fname = "testfile"
    file_path = data_dir / f"{fname}.json"
    content = {"a": 1}
    with open(file_path, "w") as f:
        json.dump(content, f)

    monkeypatch.chdir(tmp_path)

    result = read_data(fname, overrrides=["nonexistent=5"])
    assert result == content
    captured = capsys.readouterr()
    assert "Warning: Override key 'nonexistent' not found in data. Skipping." in captured.out

def test_read_data_override_type_conversion(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    fname = "testfile"
    file_path = data_dir / f"{fname}.json"
    content = {"a": 1, "b": "hello"}
    with open(file_path, "w") as f:
        json.dump(content, f)

    monkeypatch.chdir(tmp_path)

    result = read_data(fname, overrrides=["a=3.14", "b='world'"])
    assert result["a"] == 3.14
    assert result["b"] == "world"

def test_read_data_file_not_found(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError):
        read_data("doesnotexist")
