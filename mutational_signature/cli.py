import runpy


def _run(module_name):
  runpy.run_module(module_name, run_name="__main__")


def adjust_counts():
  _run("mutational_signature.adjust_counts")


def annotate_context():
  _run("mutational_signature.annotate_context")


def annotate_context_summary():
  _run("mutational_signature.annotate_context_summary")


def annotate_maf():
  _run("mutational_signature.annotate_maf")


def assign_signatures():
  _run("mutational_signature.assign_signatures")


def bootstrap():
  _run("mutational_signature.bootstrap")


def bootstrap_simple():
  _run("mutational_signature.bootstrap_simple")


def combine_counts():
  _run("mutational_signature.combine_counts")


def combine_signatures():
  _run("mutational_signature.combine_signatures")


def compare_exposures():
  _run("mutational_signature.compare_exposures")


def compare_exposures_to_signatures():
  _run("mutational_signature.compare_exposures_to_signatures")


def compare_signatures():
  _run("mutational_signature.compare_signatures")


def context_best_sig():
  _run("mutational_signature.context_best_sig")


def context_likelihood():
  _run("mutational_signature.context_likelihood")


def convert():
  _run("mutational_signature.convert")


def count():
  _run("mutational_signature.count")


def count_contexts():
  _run("mutational_signature.count_contexts")


def count_indels():
  _run("mutational_signature.count_indels")


def count_maf_indels():
  _run("mutational_signature.count_maf_indels")


def decompose():
  _run("mutational_signature.decompose")


def decompose_bulk():
  _run("mutational_signature.decompose_bulk")


def decompose_likelihood():
  _run("mutational_signature.decompose_likelihood")


def exclude_contexts():
  _run("mutational_signature.exclude_contexts")


def extended_context():
  _run("mutational_signature.extended_context")


def generate():
  _run("mutational_signature.generate")


def histograms():
  _run("mutational_signature.histograms")


def linear_dependence():
  _run("mutational_signature.linear_dependence")


def plot_components():
  _run("mutational_signature.plot_components")


def plot_counts():
  _run("mutational_signature.plot_counts")


def plot_sigs():
  _run("mutational_signature.plot_sigs")


def positive_by_chance():
  _run("mutational_signature.positive_by_chance")


def primary_context():
  _run("mutational_signature.primary_context")


def reduce_similarity():
  _run("mutational_signature.reduce_similarity")


def refine_signatures():
  _run("mutational_signature.refine_signatures")


def simulate():
  _run("mutational_signature.simulate")


def stability():
  _run("mutational_signature.stability")
