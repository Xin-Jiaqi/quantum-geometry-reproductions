"""Research handoff / protocol discovery facade, not a physics solver."""
from .registry import load_registry, validate_registry, get_module, plan_module
__all__ = ("load_registry", "validate_registry", "get_module", "plan_module")
