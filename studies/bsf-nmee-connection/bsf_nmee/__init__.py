"""BSF stacking connection -> layer-resolved NMEE tensor covariance.

Research prototype. No material-level NMEE prediction implied.
"""
from .connection import Connection, connect_interface, transform_layer_responses

__all__ = ['Connection', 'connect_interface', 'transform_layer_responses']
