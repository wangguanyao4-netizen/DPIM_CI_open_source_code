from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib import colors

# Center the plate itself in the PDF canvas. Dimension annotations occupy the
# reserved left and lower margins without shifting the plate in LaTeX.
W, H = 350, 326
x0, y0 = 70, 58
side = 210
x1, y1 = x0 + side, y0 + side
cx, cy = x0 + side/2, y0 + side/2

# colors
edge = colors.HexColor('#1f2a35')
fill = colors.HexColor('#d8e1ea')
dash = colors.HexColor('#75869a')
grid = colors.HexColor('#eef3f8')
red = colors.HexColor('#a94442')

def draw_math_label_pdf(c, x, y, parts, size=11, centered=True, angle=0):
    """Draw simple math-style label with italic variables and roman units."""
    widths = []
    for text, font_name in parts:
        c.setFont(font_name, size)
        widths.append(c.stringWidth(text, font_name, size))
    total = sum(widths)
    c.saveState()
    c.translate(x, y)
    if angle:
        c.rotate(angle)
    xx = -total / 2 if centered else 0
    for (text, font_name), width in zip(parts, widths):
        c.setFont(font_name, size)
        c.drawString(xx, 0, text)
        xx += width
    c.restoreState()

def draw_pdf(out_pdf):
    c = canvas.Canvas(str(out_pdf), pagesize=(W, H))
    c.setLineJoin(1)
    # Plate fill and border
    c.setFillColor(fill)
    c.setStrokeColor(edge)
    c.setLineWidth(1.35)
    c.rect(x0, y0, side, side, stroke=1, fill=1)
    # Subtle FEM grid
    c.setStrokeColor(grid)
    c.setLineWidth(0.25)
    n = 8
    for i in range(1, n):
        x = x0 + side*i/n
        y = y0 + side*i/n
        c.line(x, y0, x, y1)
        c.line(x0, y, x1, y)
    # Dashed simply supported boundary around four sides
    c.setStrokeColor(dash)
    c.setLineWidth(1.0)
    c.setDash(5, 4)
    off = 6
    c.rect(x0-off, y0-off, side+2*off, side+2*off, stroke=1, fill=0)
    c.setDash()
    # Center response point only; no blue load arrow
    c.setFillColor(red)
    c.circle(cx, cy, 4.0, stroke=0, fill=1)
    c.setFillColor(colors.black)
    c.setFont('Times-Italic', 10)
    c.drawCentredString(cx-2, cy-18, 'w')
    c.setFont('Times-Italic', 6.5)
    c.drawString(cx+1.8, cy-21, 'c')
    # Simply supported text
    c.setFont('Times-Bold', 13)
    c.drawCentredString(cx, cy+42, 'Simply supported')
    # Dimension lines
    c.setStrokeColor(edge)
    c.setLineWidth(0.8)
    # horizontal dimension a
    yd = y0 - 28
    c.line(x0, yd, x1, yd)
    for x,sgn in [(x0,1),(x1,-1)]:
        p=c.beginPath(); p.moveTo(x,yd); p.lineTo(x+sgn*6,yd+3); p.lineTo(x+sgn*6,yd-3); p.close(); c.drawPath(p,stroke=1,fill=1)
    c.line(x0, yd-5, x0, y0)
    c.line(x1, yd-5, x1, y0)
    draw_math_label_pdf(c, cx, yd-13, [
        ('a', 'Times-Italic'),
        (' = 1', 'Times-Roman'),
        (' m', 'Times-Roman'),
    ], size=11, centered=True)
    # vertical dimension b
    xd = x0 - 34
    c.line(xd, y0, xd, y1)
    for y,sgn in [(y0,1),(y1,-1)]:
        p=c.beginPath(); p.moveTo(xd,y); p.lineTo(xd+3,y+sgn*6); p.lineTo(xd-3,y+sgn*6); p.close(); c.drawPath(p,stroke=1,fill=1)
    c.line(xd-5, y0, x0, y0)
    c.line(xd-5, y1, x0, y1)
    draw_math_label_pdf(c, xd-13, cy, [
        ('b', 'Times-Italic'),
        (' = 1', 'Times-Roman'),
        (' m', 'Times-Roman'),
    ], size=11, centered=True, angle=90)
    c.showPage(); c.save()


