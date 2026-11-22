"""Core simulation engine: models, rate heterogeneity and evolution."""

from .aa_models import AAModel, create_aa_model
from .evolve import evolve_sequences
from .gamma import discrete_gamma, rndgamma
from .models import (
    AA_MODELS,
    NUC_MODELS,
    ModelRegistry,
    SubstitutionModel,
    build_model_registry,
    setup_model,
)
from .nuc_models import F84Model, GTRModel, HKYModel, JC69Model
from .rate_heterogeneity import RateInfo, assign_categories, normalise_branch_rates, setup_rates

__all__ = [
    "setup_model", "build_model_registry", "NUC_MODELS", "AA_MODELS",
    "SubstitutionModel", "ModelRegistry",
    "JC69Model", "HKYModel", "F84Model", "GTRModel",
    "AAModel", "create_aa_model",
    "rndgamma", "discrete_gamma",
    "RateInfo",
    "setup_rates", "assign_categories", "normalise_branch_rates",
    "evolve_sequences",
]
