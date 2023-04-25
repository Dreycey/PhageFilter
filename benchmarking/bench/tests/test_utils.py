import unittest
from collections import Counter
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
from bench.utils import compute_metrics, get_filter_metric_counts, get_classification_metric_counts,\
    get_filter_metrics, get_classification_metrics, get_readcount_metrics

class TestUtils(unittest.TestCase):
    def test_compute_metrics(self):
        # negative counts are not allowed
        with self.assertRaises(AssertionError):
            compute_metrics(-1, 53, 26)
        with self.assertRaises(AssertionError):
            compute_metrics(0, -1, 26)
        with self.assertRaises(AssertionError):
            compute_metrics(0, 53, -63)

        zero_metrics = compute_metrics(0, 0, 0)
        self.assertEqual(0, zero_metrics['recall'])
        self.assertEqual(0, zero_metrics['precision'])

        zero_tp = compute_metrics(0, 53, 72)
        self.assertEqual(0, zero_tp['recall'])
        self.assertEqual(0, zero_tp['precision'])

        example_tp = compute_metrics(90, 60, 30)
        self.assertEqual(0.75, example_tp['recall'])
        self.assertEqual(0.60, example_tp['precision'])

    def test_filter_metric_counts(self):
        """ Tests True Positive, False Positive, and False Negative counts for filtering reads. """
        true_map = {
            't1': 5729,
            't2': 6233,
            't3': 5720,
            't4': 682,
        }
        out_map = Counter({
            't1': true_map['t1'],
            't2': true_map['t2'] - 500,
            't3': true_map['t3'] + 500,
            'xyz': 5262,
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

    def test_classification_metric_counts(self):
        """ Tests True Positive, False Positive, and False Negative counts for classifying reads. """
        true_map = {
            't1': 5000,
            't2': 5000,
            't3': 5000,
            't4': 5000,
        }
        out_map = Counter({
            't1': 1,
            't2': true_map['t2'],
            't3': true_map['t3'] * 5,
            'xyz': true_map['t4'],
            'abc': 537,
        })
        metric_counts = get_classification_metric_counts(true_map, out_map)
        # t1, t2, and t3 are correctly classified. The number of reads doesn't matter currently.
        expected_tp = len(['t1', 't2', 't3'])
        expected_fp = len(['xyz', 'abc'])
        expected_fn = len(['t4'])
        self.assertEqual(expected_tp, metric_counts['TP'], "Mismatched True Positive counts")
        self.assertEqual(expected_fp, metric_counts['FP'], "Mismatched False Positive counts")
        self.assertEqual(expected_fn, metric_counts['FN'], "Mismatched False Negative counts")

    def test_filter_metrics(self):
        # The maps for this test are crafted so that the out_map perfectly counts 't1', undercounts 't2',
        # overcounts 't3', completely misses 't4', and counts the non-existent 'xyz'.
        true_map = {
            't1': 100,
            't2': 90,
            't3': 80,
            't4': 250,
        }
        out_map = Counter({
            't1': true_map['t1'],
            't2': true_map['t2'] - 10,
            't3': true_map['t3'] + 20,
            'xyz': 760,
        })
        # True Positives: 260
        # False Positives: 780
        # False Negatives: 260
        recall, precision = get_filter_metrics(true_map, out_map)
        self.assertEqual(0.5, recall)
        self.assertEqual(0.25, precision)

    def test_classification_metrics(self):
        true_map = {
            't1': 5000,
            't2': 5000,
            't3': 5000,
            't4': 5000,
        }
        out_map = Counter({
            't1': 1,
            't2': true_map['t2'],
            't3': true_map['t3'] * 5,
            'xyz': true_map['t4'],
            'abc': 537,
        })
        # True Positives: 3
        # False Positives: 2
        # False Negatives: 1
        recall, precision = get_classification_metrics(true_map, out_map)
        self.assertEqual(3/4, recall)
        self.assertEqual(3/5, precision)

    def test_readcount_metrics(self):
        true_map = {
            't1': 5729,
            't2': 6233,
            't3': 5720,
            't4': 682,
        }
        out_map = Counter({
            't1': true_map['t1'],
            't2': true_map['t2'] - 500,
            't3': true_map['t3'] + 500,
            'xyz': 5262,
        })
        # The expected differences should all be zero or positive
        expected_differences = [
            true_map['t1'] - out_map['t1'],     # Should be 0
            true_map['t2'] - out_map['t2'],     # Should be 500
            out_map['t3'] - true_map['t3'],     # Should be 500
        ]

        diffs = get_readcount_metrics(true_map, out_map)
        self.assertEqual(expected_differences, diffs)


if __name__ == '__main__':
    unittest.main()
