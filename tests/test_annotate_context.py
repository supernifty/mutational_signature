from collections import namedtuple

import mutational_signature.annotate_context as annotate_context


Variant = namedtuple("Variant", "CHROM POS REF ALT")


def test_pyrimidine_doublet_is_unchanged():
  last = Variant("1", 100, "A", ("C",))
  variant = Variant("1", 101, "C", ("A",))

  assert annotate_context.doublet_context(variant, {}, last) == "AC>CA"


def test_purine_doublet_is_reversed_and_normalized():
  last = Variant("1", 100, "A", ("C",))
  variant = Variant("1", 101, "G", ("T",))

  assert annotate_context.doublet_context(variant, {}, last) == "CT>AG"


def test_non_adjacent_variant_returns_none():
  last = Variant("1", 100, "A", ("C",))
  variant = Variant("1", 200, "G", ("T",))

  assert annotate_context.doublet_context(variant, {}, last) is None


def test_other_chromosome_returns_none():
  last = Variant("1", 100, "A", ("C",))
  variant = Variant("2", 101, "G", ("T",))

  assert annotate_context.doublet_context(variant, {}, last) is None
