"""Keep process routing and failure propagation honest at the reuse boundary."""
import subprocess
import sys
import unittest
from unittest.mock import patch
from scripts import core

class CoreEntry(unittest.TestCase):
    def test_all_registered_entries_resolve(self):
        self.assertEqual([r['id'] for r in core.capabilities()],
                         ['nmee','spin','magnetotransport','curvature-flux','igmc',
                          'altermagnet-bpve','imhe','deep-dbc','optical-design'])

    def test_routes_to_selected_module_with_current_interpreter(self):
        for row in core.capabilities():
            args,cwd,env=core.command(row,'run')
            self.assertEqual(args,[sys.executable,row['example']])
            self.assertEqual(cwd,core.ROOT/row['directory'])
            self.assertEqual(env['PYTHONPATH'],str(cwd/row['pythonpath']))

    def test_failed_solver_does_not_report_ready_result(self):
        failed=subprocess.CompletedProcess([],1,stdout='{}',stderr='solver rejected protocol')
        with patch.object(core.subprocess,'run',return_value=failed):
            with self.assertRaisesRegex(RuntimeError,'solver rejected protocol'):
                core.run_example(core.capabilities()[0])

    def test_non_object_recipe_is_rejected(self):
        bad=subprocess.CompletedProcess([],0,stdout='[]',stderr='')
        with patch.object(core.subprocess,'run',return_value=bad):
            with self.assertRaisesRegex(ValueError,'nonempty JSON object'):
                core.run_example(core.capabilities()[0])
