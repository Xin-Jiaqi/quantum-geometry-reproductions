"""Interface v1 DESIGN TYPES, not an implemented transport dispatcher.

Future adapters will implement the protocols below. This module intentionally
does not call/replace any frozen paper code. A result cannot be marked
'paper full pass' by passing this dataclass validation.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Protocol

MatrixAtK = Callable[[tuple[float, ...]], Any]

class BlochProvider(Protocol):
    """Same-orbital-gauge matrix providers; k coordinate convention explicit."""
    basis_id: str
    k_dimension: int
    k_unit: str
    matrix_unit: str
    def hamiltonian(self, k: tuple[float, ...]) -> Any: ...
    def derivative(self, axis: int, k: tuple[float, ...]) -> Any: ...
    def operator(self, name: str, k: tuple[float, ...]) -> Any: ...

@dataclass(frozen=True)
class ResponseRequestV1:
    method_id: str
    protocol_id: str
    observable_id: str
    observable_basis_id: str
    hamiltonian_basis_id: str
    occupation: str
    temperature_value: float
    temperature_unit: str
    k_convention: str
    k_unit: str
    energy_unit: str
    eta_value: float | None = None
    eta_kind: str | None = None
    frequency_value: float | None = None
    frequency_unit: str | None = None
    degeneracy_policy: str = "reject_unhandled"
    units_target: str = "model"
    source_revision: str = "UNPINNED"
    extra: Mapping[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.method_id or not self.protocol_id or not self.observable_id:
            raise ValueError("method, protocol and observable required")
        if not self.observable_basis_id or self.observable_basis_id != self.hamiltonian_basis_id:
            raise ValueError("observable and Hamiltonian MUST be in same declared gauge")
        if self.k_convention not in {"cartesian", "fractional_reciprocal"}:
            raise ValueError("coordinate convention required")
        if self.temperature_value < 0 or not self.temperature_unit:
            raise ValueError("temperature units and nonnegative T required")
        if self.eta_value is not None and self.eta_value <= 0:
            raise ValueError("eta must be > 0")
        if self.eta_kind not in (None,"adiabatic_regulator","microscopic_collision"):
            raise ValueError("eta protocol not specified")
        if self.eta_value is not None and self.eta_kind is None:
            raise ValueError("eta must identify its physical meaning")
        if self.frequency_value is not None and not self.frequency_unit:
            raise ValueError("frequency unit required")
        if self.degeneracy_policy not in {"projectors","rank_one_verified","reject_unhandled"}:
            raise ValueError("unknown degeneracy policy")
        if self.units_target == "SI" and self.energy_unit == "unspecified":
            raise ValueError("cannot claim SI from unspecified energy unit")
        if any(token in self.protocol_id for token in ("intrinsic", "tau0")) and self.eta_kind == "adiabatic_regulator":
            raise ValueError("adiabatic constant is not validated physical intrinsic DC")
        if self.source_revision == "UNPINNED":
            raise ValueError("pin the paper module source commit for a new study")

@dataclass(frozen=True)
class ResponseResultV1:
    """Proposed result provenance, not evidence that physical numbers exist."""
    protocol_id: str
    scientific_status: str
    source_commit: str
    basis_id: str
    model_units: str
    response: Any
    convergence: Mapping[str, Any]
    diagnostics: Mapping[str, Any]
    caveats: tuple[str,...]

def evaluate(_: BlochProvider, __: ResponseRequestV1) -> ResponseResultV1:
    raise NotImplementedError(
        "No unified numerical dispatcher yet: explicitly select and validate "
        "a paper-specific adapter. Consult docs/PHYSICS_INTERFACE_SPEC.md."
    )
