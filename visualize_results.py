#!/usr/bin/env python3
"""
Minimal visualize_results.py (rolled back to a simple overview and per-row images).

Reads `slot_allocations.csv` produced by `verify` and draws:
 - one PNG per row: output/row_<r>.png
 - an overview stacking rows: output/slots_overview.png

This file intentionally keeps behavior small and dependency-light (only matplotlib required).
"""
import os
import csv
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib as mpl

# Prefer macOS CJK-capable fonts for Chinese labels; fall back to DejaVu Sans
mpl.rcParams['font.sans-serif'] = [
    'PingFang SC', 'Heiti TC', 'Hiragino Sans GB', 'STHeiti', 'Arial Unicode MS', 'DejaVu Sans'
]
mpl.rcParams['axes.unicode_minus'] = False


def find_slot_alloc_csv(base='output'):
    # prefer direct file under base, otherwise search output/output_* subdirs
    direct = os.path.join(base, 'slot_allocations.csv')
    if os.path.exists(direct):
        return direct
    if not os.path.exists(base):
        return None
    for name in os.listdir(base):
        path = os.path.join(base, name)
        if os.path.isdir(path) and name.startswith('output'):
            candidate = os.path.join(path, 'slot_allocations.csv')
            if os.path.exists(candidate):
                return candidate
    return None


def load_slots(path):
    # returns list of records (dict), and inferred numRows, numSlotsPerRow
    recs = []
    max_row = -1
    max_slot = -1
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            try:
                ship = int(r['ship'])
                k = int(r.get('k', -1)) if r.get('k', '') != '' else -1
                row = int(r['row'])
                slot = int(r['slot'])
            except Exception:
                continue
            recs.append({'ship': ship, 'k': k, 'row': row, 'slot': slot})
            if row > max_row:
                max_row = row
            if slot > max_slot:
                max_slot = slot
    numRows = max_row + 1 if max_row >= 0 else 0
    numSlotsPerRow = max_slot + 1 if max_slot >= 0 else 0
    return recs, numRows, numSlotsPerRow


def build_grid(recs, numRows, numSlotsPerRow):
    # grid stores either -1 (empty) or a tuple (ship, k)
    grid = [[-1 for _ in range(numSlotsPerRow)] for _ in range(numRows)]
    ships = set()
    for r in recs:
        s = r['ship']
        k = r['k']
        row = r['row']
        slot = r['slot']
        ships.add(s)
        if 0 <= row < numRows and 0 <= slot < numSlotsPerRow:
            grid[row][slot] = (s, k)
    return grid, sorted(list(ships))


def ship_colors(ships):
    # generate distinct colors for ships
    cmap = plt.get_cmap('tab20')
    colors = {}
    for i, s in enumerate(ships):
        colors[s] = cmap(i % 20)
    return colors


def _text_color_for_bg(color):
    # Choose black/white text for contrast given an RGB(A) color tuple (0-1)
    try:
        r, g, b = color[0], color[1], color[2]
        lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
        return 'white' if lum < 0.6 else 'black'
    except Exception:
        return 'black'


def draw_row(grid, r, ship_to_color, outdir):
    numSlots = len(grid[r])
    fig, ax = plt.subplots(figsize=(max(6, numSlots / 4), 1.2))
    for v in range(numSlots):
        cell = grid[r][v]
        if cell == -1:
            color = (0.95, 0.95, 0.95)
            rect = plt.Rectangle((v, 0), 1, 1, facecolor=color, edgecolor='k', linewidth=0.3)
            ax.add_patch(rect)
        else:
            s, k = cell
            color = ship_to_color.get(s, (0.6, 0.6, 0.6))
            rect = plt.Rectangle((v, 0), 1, 1, facecolor=color, edgecolor='k', linewidth=0.3)
            ax.add_patch(rect)
            # label with ship and compartment
            # choose a readable font size according to number of slots
            fontsize = max(6, int(12 - numSlots / 4))
            # center text
            ax.text(v + 0.5, 0.5, f'S{s}K{k}', ha='center', va='center', fontsize=fontsize, color=_text_color_for_bg(color))
    ax.set_xlim(0, numSlots)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title(f'Row {r}')
    plt.tight_layout()
    fname = os.path.join(outdir, f'row_{r:02d}.png')
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    return fname


def draw_overview(grid, ship_to_color, outdir):
    numRows = len(grid)
    numSlots = len(grid[0]) if numRows > 0 else 0
    fig_h = max(6, numSlots / 4)
    fig, axes = plt.subplots(numRows, 1, figsize=(fig_h, max(2, numRows * 0.6)))
    if numRows == 1:
        axes = [axes]
    for r, ax in enumerate(axes):
        for v in range(numSlots):
            cell = grid[r][v]
            if cell == -1:
                color = (0.95, 0.95, 0.95)
                rect = plt.Rectangle((v, 0), 1, 1, facecolor=color, edgecolor='k', linewidth=0.3)
                ax.add_patch(rect)
            else:
                s, k = cell
                color = ship_to_color.get(s, (0.6, 0.6, 0.6))
                rect = plt.Rectangle((v, 0), 1, 1, facecolor=color, edgecolor='k', linewidth=0.3)
                ax.add_patch(rect)
                fontsize = max(5, int(10 - numSlots / 6))
                ax.text(v + 0.5, 0.5, f'S{s}K{k}', ha='center', va='center', fontsize=fontsize, color=_text_color_for_bg(color))
        ax.set_xlim(0, numSlots)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel(f'R{r}', rotation=0, labelpad=20)
    plt.suptitle('Slot allocations overview (rows stacked)')
    # prepare legend patches and place legend to the right
    try:
        ships_sorted = sorted(ship_to_color.keys())
        patches = [mpatches.Patch(color=ship_to_color[s], label=f'Ship {s}') for s in ships_sorted]
        fig.subplots_adjust(right=0.78)
        fig.legend(handles=patches, bbox_to_anchor=(0.99, 0.5), loc='center left', fontsize=8)
    except Exception:
        pass
    plt.tight_layout()
    out = os.path.join(outdir, 'slots_overview.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return out


def save_legend(ship_to_color, outdir):
    # create legend patches for ships
    patches = []
    ships = sorted(ship_to_color.keys())
    for s in ships:
        patches.append(mpatches.Patch(color=ship_to_color[s], label=f'Ship {s}'))
    try:
        lg_fig = plt.figure(figsize=(2.5, max(1.5, len(ships) * 0.25)))
        lg_ax = lg_fig.add_subplot(111)
        lg_ax.axis('off')
        lg = lg_ax.legend(handles=patches, loc='center')
        for text in lg.get_texts():
            text.set_fontsize(8)
        lg_fig.savefig(os.path.join(outdir, 'legend.png'), bbox_inches='tight', dpi=150)
        plt.close(lg_fig)
    except Exception:
        pass


def find_csv_under_output(name, base='output'):
    direct = os.path.join(base, name)
    if os.path.exists(direct):
        return direct
    if not os.path.exists(base):
        return None
    for nm in os.listdir(base):
        path = os.path.join(base, nm)
        if os.path.isdir(path) and nm.startswith('output'):
            cand = os.path.join(path, name)
            if os.path.exists(cand):
                return cand
    return None


def plot_berth_assignment(outdir, ship_to_color):
    berth_csv = find_csv_under_output('berth_assignment.csv', outdir)
    if not berth_csv:
        return None
    berth_rows = []
    with open(berth_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            try:
                ship = int(rec.get('ship', '').strip())
                berth_start = float(rec.get('berth_start', '').strip())
                berth_end = float(rec.get('berth_end', '').strip())
            except Exception:
                continue
            berth_rows.append((ship, berth_start, berth_end))
    if not berth_rows:
        return None
    # sort by ship id for stable ordering
    berth_rows.sort(key=lambda x: x[0])
    ships = [s for s, _, _ in berth_rows]
    n = len(ships)
    fig, ax = plt.subplots(figsize=(8, max(3, n * 0.6)))
    # draw bars from berth_start to berth_end
    bar_height = 0.6
    for i, (s, b_start, b_end) in enumerate(berth_rows):
        color = ship_to_color.get(s, (0.6, 0.6, 0.6))
        bar_width = b_end - b_start
        rect = plt.Rectangle((b_start, i - bar_height / 2), bar_width, bar_height, facecolor=color, edgecolor='k')
        ax.add_patch(rect)
        ax.text(b_start + bar_width / 2, i, f'Ship {s}', ha='center', va='center', color=_text_color_for_bg(color), fontsize=9)
    # y axis setup with padding to fully show top/bottom labels
    ax.set_yticks(list(range(n)))
    ax.set_yticklabels([f'Ship {s}' for s in ships])
    ax.set_ylim(-0.6, n - 0.4)
    # x axis为泊位起始和终止坐标
    all_starts = [b for _, b, _ in berth_rows]
    all_ends = [e for _, _, e in berth_rows]
    # 横坐标范围固定为0-1000，每50显示一个刻度
    ax.set_xlim(0, 1000)
    ax.set_xticks(range(0, 1001, 50))
    ax.set_xlabel('Berth position')
    ax.set_title('Berth assignment (起止位置)')
    fig.subplots_adjust(left=0.18, right=0.96, top=0.88, bottom=0.18)
    out = os.path.join(outdir, 'berth_assignment.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return out


def plot_start_times(outdir, ship_to_color):
    e_s_csv = find_csv_under_output('e_s.csv', outdir)
    e_sk_csv = find_csv_under_output('e_sk.csv', outdir)
    if not e_s_csv or not e_sk_csv:
        return None
    e_s_rows = []
    with open(e_s_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            try:
                ship = int(rec.get('ship', '').strip())
                e_s = float(rec.get('e_s', '').strip())
            except Exception:
                continue
            e_s_rows.append({'ship': ship, 'e_s': e_s})
    e_sk_rows = []
    with open(e_sk_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            try:
                ship = int(rec.get('ship', '').strip())
                k = int(rec.get('k', '').strip()) if rec.get('k', '').strip() != '' else None
                e_sk = float(rec.get('e_sk', '').strip())
            except Exception:
                continue
            e_sk_rows.append({'ship': ship, 'k': k, 'e_sk': e_sk})
    ships_all = sorted({r['ship'] for r in e_s_rows})
    ship_index = {s: i for i, s in enumerate(ships_all)}
    fig, ax = plt.subplots(figsize=(10, max(3, len(ships_all) * 0.6)))
    for r in e_s_rows:
        s = r['ship']
        t = r['e_s']
        # 船级别的开始时间点，放在该船的固定y位置，保证横向对齐
        ax.scatter(t, ship_index[s], marker='D', s=80, color=ship_to_color.get(s, (0.2, 0.2, 0.8)))
    # group e_sk by ship
    es_by_ship = {}
    for r in e_sk_rows:
        es_by_ship.setdefault(r['ship'], []).append(r)
    for s, recs in es_by_ship.items():
        recs_sorted = sorted(recs, key=lambda x: x['e_sk'])
        for idx, rec in enumerate(recs_sorted):
            t = rec['e_sk']
            k = rec['k']
            y = ship_index[s]
            # 货舱开始时间与该船同一条水平线上，点置于 (t, y)
            ax.scatter(t, y, marker='o', s=30, color='black')
            if k is not None:
                # 将标注放在点的右侧，使用像素偏移，避免随坐标缩放
                ax.annotate(f'k{k}', xy=(t, y), xytext=(4, 0), textcoords='offset points',
                            ha='left', va='center', fontsize=7)
    ax.set_yticks(list(ship_index.values()))
    ax.set_yticklabels([f'Ship {s}' for s in ships_all])
    ax.set_xlabel('Time (hours)')
    ax.set_title('Start times: e_s (diamond) and e_sk (circle)')
    # 加细横向参考线以增强对齐感知
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    # Place legend outside the axes to avoid covering points
    try:
        handles = [
            Line2D([0], [0], marker='D', color='none', markerfacecolor=(0.2, 0.2, 0.8), markersize=8, label='e_s (ship start)'),
            Line2D([0], [0], marker='o', color='black', markersize=6, linestyle='None', label='e_sk (compartment)')
        ]
        fig.subplots_adjust(right=0.78)
        fig.legend(handles=handles, bbox_to_anchor=(0.99, 0.5), loc='center left', fontsize=8)
    except Exception:
        pass
    timeline_out = os.path.join(outdir, 'start_times.png')
    plt.tight_layout()
    fig.savefig(timeline_out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return timeline_out


def main():
    # 自动在 output/ 或 output/output_* 下寻找最近的 slot_allocations.csv
    csv_path = find_slot_alloc_csv('output/output_2')
    if csv_path is None:
        print('Could not find slot_allocations.csv under output/ or output/output_*/')
        return
    recs, numRows, numSlotsPerRow = load_slots(csv_path)
    if numRows == 0 or numSlotsPerRow == 0:
        print('No slot data found in', csv_path)
        return
    grid, ships = build_grid(recs, numRows, numSlotsPerRow)
    ship_to_color = ship_colors(ships)
    # Save images next to the detected CSV (e.g., output/output_1)
    outdir = os.path.dirname(csv_path) if csv_path else 'output'
    os.makedirs(outdir, exist_ok=True)
    imgs = []
    for r in range(numRows):
        imgs.append(draw_row(grid, r, ship_to_color, outdir))
    overview = draw_overview(grid, ship_to_color, outdir)
    print('Generated', overview)
    for p in imgs[:5]:
        print('Row image:', p)
    # save legend
    lg = save_legend(ship_to_color, outdir)
    if lg is not None:
        print('Generated', os.path.join(outdir, 'legend.png'))
    # berth assignment
    ba = plot_berth_assignment(outdir, ship_to_color)
    if ba:
        print('Generated', ba)
    # start times
    st = plot_start_times(outdir, ship_to_color)
    if st:
        print('Generated', st)


if __name__ == '__main__':
    main()
