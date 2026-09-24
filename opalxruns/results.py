"""The result table the scripted studies print into results.txt."""

from __future__ import annotations


class Results:
    """Rows of (test, quantity, measured, expected, tol) printed with PASS/FAIL.
    A row with tol None is a diagnostic and never fails the suite.

    The keyword arguments only change how print_table() lays the table out, so
    each study keeps the results.txt it has always written:
      first, first_width, first_align   header, width and alignment of column 1
      name_width                        width of the quantity column
      expected                          header of the expected-value column
      rule                              length of the '=' lines around a title
    """

    def __init__(self, first="test", first_width=4, first_align=">", name_width=34,
                 expected="expected", rule=104):
        self.rows = []
        self.first = first
        self.first_width = first_width
        self.first_align = first_align
        self.name_width = name_width
        self.expected = expected
        self.rule = rule

    def check(self, test, name, measured, expected, tol, rel=False, note=""):
        if tol is None:
            ok = None
        elif rel:
            ok = bool(abs(measured - expected) / max(abs(expected), 1e-30) <= tol)
        else:
            ok = bool(abs(measured - expected) <= tol)  # bool(): numpy bools fail `is False`
        self.rows.append(dict(test=test, name=name, measured=measured, expected=expected,
                              tol=tol, rel=rel, ok=ok, note=note))
        return ok

    def at_least(self, test, name, measured, minimum, note=""):
        """A check that something is big enough rather than small enough, for
        the cases where the point is that two results must NOT agree."""
        self.rows.append(dict(test=test, name=name, measured=measured,
                              expected=minimum, tol=None, rel=False,
                              ok=bool(measured >= minimum),
                              note=note or f"must be at least {minimum:g}"))
        return measured >= minimum

    def note(self, test, name, text):
        self.rows.append(dict(test=test, name=name, measured=None, expected=None,
                              tol=None, rel=False, ok=None, note=text))

    @property
    def failed(self) -> int:
        return sum(1 for r in self.rows if r["ok"] is False)

    def print_table(self, title=""):
        w, a, nw = self.first_width, self.first_align, self.name_width
        if title:
            print(f"\n{'=' * self.rule}\n{title}\n{'=' * self.rule}")
        hdr = (f"{self.first:{w}s} {'quantity':{nw}s} {'measured':>13s} {self.expected:>13s} "
               f"{'|diff|':>10s} {'tol':>9s}  result")
        print(hdr)
        print("-" * len(hdr))
        for r in self.rows:
            if r["measured"] is None:
                print(f"{r['test']:{a}{w}} {r['name']:{nw}s} {r['note']}")
                continue
            diff = abs(r["measured"] - r["expected"])
            res = "diag" if r["ok"] is None else ("PASS" if r["ok"] else "**FAIL**")
            tolstr = "-" if r["tol"] is None else f"{r['tol']:.1e}{'r' if r['rel'] else ''}"
            note = f"  {r['note']}" if r["note"] else ""
            print(f"{r['test']:{a}{w}} {r['name']:{nw}s} {r['measured']:>13.6g} "
                  f"{r['expected']:>13.6g} {diff:>10.3g} {tolstr:>9s}  {res}{note}")
