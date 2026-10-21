"""python -m integration.physics_interface.cli list|inspect|plan|validate.

Read-only protocol discovery. Never spawns subprocesses.
"""
from __future__ import annotations
import argparse
import json
from .registry import (InvalidRegistry, get_module, load_registry,
                       plan_module, validate_registry)

def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description="Inspect versioned paper method contracts")
    cmds=parser.add_subparsers(dest="command",required=True)
    cmds.add_parser("list")
    cmds.add_parser("validate")
    item=cmds.add_parser("inspect");item.add_argument("module")
    plan=cmds.add_parser("plan");plan.add_argument("module");plan.add_argument("protocol")
    args=parser.parse_args(argv)
    try:
        reg=load_registry()
        if args.command=="list":
            out=[{"id":x["id"],"paper":x["paper"],"status":x["science_status"],
                  "protocol_ids":[p["id"] for p in x["protocols"]]}
                 for x in reg["modules"]]
        elif args.command=="validate":
            out=validate_registry(reg,check_files=True)
        elif args.command=="inspect":
            out=get_module(args.module,reg)
        else:
            out=plan_module(args.module,args.protocol,reg)
    except (InvalidRegistry, KeyError, OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
