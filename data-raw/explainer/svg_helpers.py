"""Small SVG helpers shared by the explainer panels."""
import math, random

def esc(s): return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
def T(x, y, s, cls="tm", a="start", extra=""):
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{a}" {extra}>{esc(s)}</text>'
def R(x, y, w, h, cls, rx=3, extra=""):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" class="{cls}" {extra}/>'
def C(x, y, r, cls, extra=""):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" class="{cls}" {extra}/>'
def L(x1, y1, x2, y2, cls, extra=""):
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="{cls}" {extra}/>'
def P(d, cls, extra=""):
    return f'<path d="{d}" class="{cls}" {extra}/>'
def arrow(x1, y1, x2, y2, cls="ar"):
    ang = math.atan2(y2 - y1, x2 - x1); h = 6
    a1 = (x2 - h * math.cos(ang - .45), y2 - h * math.sin(ang - .45))
    a2 = (x2 - h * math.cos(ang + .45), y2 - h * math.sin(ang + .45))
    return (L(x1, y1, x2, y2, cls) +
            P(f"M{a1[0]:.1f} {a1[1]:.1f} L{x2:.1f} {y2:.1f} L{a2[0]:.1f} {a2[1]:.1f}", cls))
def badge(x, y, n, optional=False, color="1"):
    """numbered step badge centred at (x, y)"""
    if optional:
        return C(x, y, 8, f"fsurf st{color} w15") + T(x, y + 3.5, str(n), f"tny f{color}t", "middle")
    return C(x, y, 8, f"f{color}") + T(x, y + 3.5, str(n), "tny finv", "middle")

def svg(w, h, body, label):
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{esc(label)}">{body}</svg>'
def wave(x, y, w, amp, n=14, cls="wv"):
    pts = []
    rnd = random.Random(int(x * 7 + y))
    for i in range(n + 1):
        xx = x + w * i / n
        a = amp * (0.25 + 0.75 * math.sin(math.pi * i / n)) * (1 if i % 2 else -1) * (0.6 + 0.4 * rnd.random())
        pts.append(f"{xx:.1f} {y + a:.1f}")
    return P("M" + " L".join(pts), cls)


NOTE = 'font-size="10.5" opacity=".85"'
MONO = 'style="font-size:8.4px"'
def note(b, x, y, s): b.append(T(x, y, s, "ts", "start", NOTE))
