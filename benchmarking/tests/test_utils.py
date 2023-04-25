import unittest

from collections import Counter

from bench.utils import get_filter_metric_counts


class TestUtils(unittest.TestCase):
    def test_filter_metrics(self):
        """ Tests True Positive, False Positive, and False Negative counts for filtering reads. """
        true_map = {
            't1': 5729,
            't2': 6233,
            't3': 5720,
            't4': 682
        }
        out_map = Counter({
            't1': true_map['t1'],
            't2': true_map['t2'] - 500,
            't3': true_map['t3'] + 500,
            'xyz': 5262
        })
        metric_counts = get_filter_metric_counts(true_map, out_map)
        # The maps for this test are crafted so that the out_map perfectly counts 't1', undercounts 't2',
        # overcounts 't3', completely misses 't4', and counts the non-existent 'xyz'.
        expected_tp = true_map['t1'] + out_map['t2'] + true_map['t3']
        expected_fp = (out_map['t3'] - true_map['t3']) + out_map['xyz']
        expected_fn = (true_map['t2'] - out_map['t2']) + true_map['t4']
        self.assertEqual(expected_tp, metric_counts['TP'], "Mismatched True Positive counts")
        self.assertEqual(expected_fp, metric_counts['FP'], "Mismatched False Positive counts")
        self.assertEqual(expected_fn, metric_counts['FN'], "Mismatched False Negative counts")


if __name__ == '__main__':
    unittest.main()
