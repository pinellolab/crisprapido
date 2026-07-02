#!/usr/bin/env python3
"""
Generate test cases by running the Python reference CFD implementation.
This script generates test cases with expected values from the reference implementation.
"""

import sys
import os
import json

# Add the reference Python implementation to path
script_dir = os.path.dirname(os.path.abspath(__file__))
reference_dir = os.path.join(script_dir, '../../reference_cfd/python')
sys.path.insert(0, reference_dir)
os.chdir(reference_dir)

from cfd_score_calculator import calculate_cfd

def generate_test_cases():
    """Generate comprehensive test cases with expected values from Python reference"""
    test_cases = []

    base_seq = "GAAACAGTCGATTTTATCAC"  # 20nt test sequence
    bases = ['A', 'C', 'G', 'T']
    pams = ['GG', 'AG', 'CG', 'TG', 'GA', 'AA', 'CA', 'TA']

    # 1. Perfect match with various PAMs
    for pam in ['GG', 'AG', 'CG', 'TG', 'GA', 'AA', 'CA', 'TA', 'GT', 'AT', 'CT', 'TT',
                'GC', 'AC', 'CC', 'TC']:
        try:
            score = calculate_cfd(base_seq, base_seq, pam)
            test_cases.append({
                'spacer': base_seq,
                'protospacer': base_seq,
                'pam': pam,
                'expected_score': score,
                'description': f'Perfect match with PAM {pam}'
            })
        except Exception as e:
            test_cases.append({
                'spacer': base_seq,
                'protospacer': base_seq,
                'pam': pam,
                'expected_error': str(e),
                'description': f'Perfect match with PAM {pam}'
            })

    # 2. Single mismatches at each position (1-20)
    for pos in range(20):
        original = base_seq[pos]
        for new_base in bases:
            if new_base == original:
                continue
            proto = base_seq[:pos] + new_base + base_seq[pos+1:]
            try:
                score = calculate_cfd(base_seq, proto, 'GG')
                test_cases.append({
                    'spacer': base_seq,
                    'protospacer': proto,
                    'pam': 'GG',
                    'expected_score': score,
                    'description': f'Mismatch at position {pos+1}: {original}->{new_base}'
                })
            except Exception as e:
                test_cases.append({
                    'spacer': base_seq,
                    'protospacer': proto,
                    'pam': 'GG',
                    'expected_error': str(e),
                    'description': f'Mismatch at position {pos+1}: {original}->{new_base}'
                })

    # 3. Gap in spacer at each position (2-20, position 1 is free)
    for pos in range(20):
        spacer_with_gap = base_seq[:pos] + '-' + base_seq[pos+1:]
        try:
            score = calculate_cfd(spacer_with_gap, base_seq, 'GG')
            test_cases.append({
                'spacer': spacer_with_gap,
                'protospacer': base_seq,
                'pam': 'GG',
                'expected_score': score,
                'description': f'Gap in spacer at position {pos+1}'
            })
        except Exception as e:
            test_cases.append({
                'spacer': spacer_with_gap,
                'protospacer': base_seq,
                'pam': 'GG',
                'expected_error': str(e),
                'description': f'Gap in spacer at position {pos+1}'
            })

    # 4. Gap in protospacer at each position (2-20, position 1 is free)
    for pos in range(20):
        proto_with_gap = base_seq[:pos] + '-' + base_seq[pos+1:]
        try:
            score = calculate_cfd(base_seq, proto_with_gap, 'GG')
            test_cases.append({
                'spacer': base_seq,
                'protospacer': proto_with_gap,
                'pam': 'GG',
                'expected_score': score,
                'description': f'Gap in protospacer at position {pos+1}'
            })
        except Exception as e:
            test_cases.append({
                'spacer': base_seq,
                'protospacer': proto_with_gap,
                'pam': 'GG',
                'expected_error': str(e),
                'description': f'Gap in protospacer at position {pos+1}'
            })

    # 5. Multiple mismatches (2-3 at various positions)
    multi_mismatch_cases = [
        (1, 20),   # PAM-distal and PAM-proximal
        (1, 8),    # PAM-distal and middle
        (8, 18),   # Two middle positions
        (5, 10, 15),  # Three positions
        (1, 10, 20),  # Three spanning positions
    ]

    for positions in multi_mismatch_cases:
        proto = list(base_seq)
        desc_parts = []
        for pos in positions:
            original = proto[pos-1]
            # Change to next base cyclically
            new_base = bases[(bases.index(original) + 1) % 4]
            proto[pos-1] = new_base
            desc_parts.append(f'pos{pos}:{original}->{new_base}')
        proto = ''.join(proto)
        try:
            score = calculate_cfd(base_seq, proto, 'GG')
            test_cases.append({
                'spacer': base_seq,
                'protospacer': proto,
                'pam': 'GG',
                'expected_score': score,
                'description': f'Multiple mismatches: {", ".join(desc_parts)}'
            })
        except Exception as e:
            test_cases.append({
                'spacer': base_seq,
                'protospacer': proto,
                'pam': 'GG',
                'expected_error': str(e),
                'description': f'Multiple mismatches: {", ".join(desc_parts)}'
            })

    # 6. Mismatch + non-canonical PAM
    for pam in ['AG', 'CG', 'TG']:
        # Single mismatch at position 10
        proto = base_seq[:9] + ('C' if base_seq[9] != 'C' else 'A') + base_seq[10:]
        try:
            score = calculate_cfd(base_seq, proto, pam)
            test_cases.append({
                'spacer': base_seq,
                'protospacer': proto,
                'pam': pam,
                'expected_score': score,
                'description': f'Mismatch at pos 10 with PAM {pam}'
            })
        except Exception as e:
            test_cases.append({
                'spacer': base_seq,
                'protospacer': proto,
                'pam': pam,
                'expected_error': str(e),
                'description': f'Mismatch at pos 10 with PAM {pam}'
            })

    # 7. Case sensitivity tests
    case_tests = [
        (base_seq.lower(), base_seq, 'GG', 'lowercase spacer'),
        (base_seq, base_seq.lower(), 'GG', 'lowercase protospacer'),
        (base_seq.lower(), base_seq.lower(), 'gg', 'all lowercase'),
    ]

    for spacer, proto, pam, desc in case_tests:
        try:
            score = calculate_cfd(spacer, proto, pam)
            test_cases.append({
                'spacer': spacer,
                'protospacer': proto,
                'pam': pam,
                'expected_score': score,
                'description': desc
            })
        except Exception as e:
            test_cases.append({
                'spacer': spacer,
                'protospacer': proto,
                'pam': pam,
                'expected_error': str(e),
                'description': desc
            })

    # 8. Real-world examples from literature
    literature_cases = [
        ("CTAACAGTTGCTTTTATCAC", "CTAACAGTTGCTTTTATCAC", "GG", "Perfect on-target"),
        ("CTAACAGTTGCTTTTATCAC", "TTAACAGTTGCTTTTATCAC", "GG", "Literature: pos1 mismatch"),
        ("CTAACAGTTGCTTTTATCAC", "CTAACAGATGCTTTTATCAC", "GG", "Literature: pos7 mismatch"),
    ]

    for spacer, proto, pam, desc in literature_cases:
        try:
            score = calculate_cfd(spacer, proto, pam)
            test_cases.append({
                'spacer': spacer,
                'protospacer': proto,
                'pam': pam,
                'expected_score': score,
                'description': desc
            })
        except Exception as e:
            test_cases.append({
                'spacer': spacer,
                'protospacer': proto,
                'pam': pam,
                'expected_error': str(e),
                'description': desc
            })

    return test_cases

def main():
    test_cases = generate_test_cases()

    # Output as JSON
    output_file = os.path.join(script_dir, 'expected_values.json')
    with open(output_file, 'w') as f:
        json.dump(test_cases, f, indent=2)

    print(f"Generated {len(test_cases)} test cases")
    print(f"Output written to: {output_file}")

    # Also print summary
    with_scores = sum(1 for tc in test_cases if 'expected_score' in tc)
    with_errors = sum(1 for tc in test_cases if 'expected_error' in tc)
    print(f"  - Cases with scores: {with_scores}")
    print(f"  - Cases with errors: {with_errors}")

if __name__ == '__main__':
    main()
