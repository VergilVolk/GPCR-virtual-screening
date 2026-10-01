"""Compatibility entry point for the corrected PACER-200 handoff audit.

The former implementation compared PACER-200 against a different 28,519-
molecule screening library and therefore could not measure pipeline continuity.
"""
from audit_candidate_handoff_v02 import main


if __name__ == "__main__":
    main()
