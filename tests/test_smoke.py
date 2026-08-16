import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


fleet = load_module("bakta_fleet", "scripts/anot00_bakta_01_fleet.py")
selector = load_module("bakta_selector", "scripts/anot00_bakta_02_selecionados.py")


class FleetSmokeTests(unittest.TestCase):
    def test_primary_genomic_fasta_is_accepted(self):
        self.assertTrue(fleet.is_primary_genomic_fna("GCF_000001_genomic.fna"))

    def test_cds_and_rna_fastas_are_rejected(self):
        self.assertFalse(
            fleet.is_primary_genomic_fna("GCF_000001_cds_from_genomic.fna")
        )
        self.assertFalse(
            fleet.is_primary_genomic_fna("GCF_000001_rna_from_genomic.fna")
        )

    def test_sample_name_is_sanitized(self):
        self.assertEqual(fleet.safe_sample_name("Sample 01!"), "Sample_01")

    def test_accession_normalization_removes_version_and_nz_prefix(self):
        self.assertEqual(fleet.normalize_accession("NZ_CP012345.2"), "CP012345")


class SelectorSmokeTests(unittest.TestCase):
    def setUp(self):
        self.params = selector.SelectionParameters(
            short_aa_threshold=90,
            weak_identity=0.90,
            weak_query_cov=0.80,
            weak_subject_cov=0.80,
            include_sorf_targets=False,
            include_pseudogenes=False,
            allow_id_mismatch=False,
            allow_sequence_mismatch=False,
            allow_hypothetical_mismatch=False,
        )

    def test_hypothetical_product_is_high_priority(self):
        inference = selector.InferenceRecord()
        result = selector.classify_target(
            "hypothetical protein", inference, self.params
        )
        self.assertIsNotNone(result)
        self.assertEqual(result[0], "high")

    def test_generic_domain_product_is_medium_priority(self):
        inference = selector.InferenceRecord()
        result = selector.classify_target(
            "DUF1234 family protein", inference, self.params
        )
        self.assertIsNotNone(result)
        self.assertEqual(result[0], "medium")

    def test_uncertain_product_is_low_priority(self):
        inference = selector.InferenceRecord()
        result = selector.classify_target(
            "putative transporter", inference, self.params
        )
        self.assertIsNotNone(result)
        self.assertEqual(result[0], "low")

    def test_weak_identity_is_low_priority(self):
        inference = selector.InferenceRecord(
            method="psc",
            accession="UniRef90_example",
            identity=0.75,
            query_cov=0.95,
            subject_cov=0.95,
        )
        result = selector.classify_target(
            "ABC transporter ATP-binding protein", inference, self.params
        )
        self.assertIsNotNone(result)
        self.assertEqual(result[0], "low")
        self.assertIn("weak_identity", result[1])

    def test_strong_informative_annotation_is_not_selected(self):
        inference = selector.InferenceRecord(
            method="psc",
            accession="UniRef90_example",
            identity=0.97,
            query_cov=0.95,
            subject_cov=0.95,
        )
        result = selector.classify_target(
            "ABC transporter ATP-binding protein", inference, self.params
        )
        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
