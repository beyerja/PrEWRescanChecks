""" Helper function related to the system that your running on.
"""

from pathlib import Path
import os

def create_dir(dir):
  """ Try to create the given directory and it's parents.
  """
  Path(dir).mkdir(parents=True, exist_ok=True)
