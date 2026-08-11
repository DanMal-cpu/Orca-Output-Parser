# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 13:28:57 2025

@author: dmalcolm
"""

import numpy as np

from matplotlib.pyplot import savefig, show, subplots, tight_layout
from os import path
from re import match
from scipy.signal import find_peaks




max_peaks_to_label = 5
show_transmittance = False
sigma_value = 25
stick_colour = 'gray'
stick_alpha = 0.7
stick_linewidth = 1.2
stick_max_relative_height = 0.5
timeout = 10


def add_spectrum_annotations(ax, peaks, colors, stick_wavelengths=None, stick_fosc=None, max_fosc=1.0, label_prefix='λ =', y_offset_base=15):
    """ Add annotations for peaks on a spectrum, automatically avoiding stick overlap """
    
    if stick_wavelengths is None or stick_fosc is None:
        for i, (wl, val) in enumerate(peaks[:5]):
            ax.plot(wl, val, 'o', color=colors[i], markersize=8, markeredgecolor='black', markeredgewidth=1)
            
            if val > 0.8:
                x_offset, y_offset, va, ha = 0, -30, 'top', 'center'
            else:
                x_offset = 10 if i % 2 == 0 else -10
                y_offset = y_offset_base if i % 3 == 0 else -y_offset_base
                va = 'bottom' if y_offset > 0 else 'top'
                ha = 'center' if val > 0.8 else ('left' if x_offset > 0 else 'right')
            
            ax.annotate(f'{label_prefix} {wl:.1f} nm', xy=(wl, val), xytext=(x_offset, y_offset), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow' if 'Abs' in label_prefix else 'lightyellow', alpha=0.9, edgecolor='black'), fontsize=14, ha=ha, va=va)
        return
    
    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    points_to_data = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (bbox.height * 72) if bbox.height > 0 else 0.0025
    
    for i, (wl, val) in enumerate(peaks[:5]):
        ax.plot(wl, val, 'o', color=colors[i], markersize=8, markeredgecolor='black', markeredgewidth=1)
        
        if val > 0.8:
            x_offset, y_offset, va, ha = 0, -30, 'top', 'center'
        else:
            x_offset = 10 if i % 2 == 0 else -10
            y_offset = y_offset_base if i % 3 == 0 else -y_offset_base
            va = 'bottom' if y_offset > 0 else 'top'
            ha = 'center' if val > 0.8 else ('left' if x_offset > 0 else 'right')
        
        move_above = False
        if y_offset < 0:
            box_bottom = val + y_offset * points_to_data
            for stick_wl, stick_f in zip(stick_wavelengths, stick_fosc):
                if abs(wl - stick_wl) < 25.0:
                    stick_height = (stick_f / max_fosc) * stick_max_relative_height * 1.0
                    if stick_height > box_bottom:
                        move_above = True
                        break
        
        if move_above:
            x_offset, y_offset, va, ha = 0, 30, 'bottom', 'center'
        
        ax.annotate(f'{label_prefix} {wl:.1f} nm', xy=(wl, val), xytext=(x_offset, y_offset), 
                   textcoords='offset points', 
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow' if 'Abs' in label_prefix else 'lightyellow', alpha=0.9, edgecolor='black'), fontsize=14, ha=ha, va=va)



def find_major_peaks(data, wavelengths, min_prominence=0.05, min_height=0.1):
    """ Find major peaks in spectral data with customizable thresholds """
    
    peaks, _ = find_peaks(data, prominence=min_prominence, height=min_height)
    
    peak_wavelengths = wavelengths[peaks]
    peak_values = data[peaks]
    
    sorted_indices = np.argsort(peak_values)[::-1]          # Sort peaks by value (descending) to get major peaks first.
    return list(zip(peak_wavelengths[sorted_indices], peak_values[sorted_indices]))



def generate_spectrum(filename, show_transmittance, sigma_value):
    """ Main function to analyze TD-DFT spectrum from a file with interactive sigma adjustment. 
    Returns png plot and dict of major peaks"""
    
    
    # ============ PARSE TD-DFT DATA FROM FILE USING GIVEN SIGMA VALUE ============
    
    tddft_peaks = parse_tddft_file(filename, sigma_value)
    
    if not tddft_peaks:
        print(f"No valid TD-DFT data found in {filename}")
        return None
    
    print(f"Parsed TD-DFT Peaks From {filename}\n")


    # ============ DEFINE THE WAVELENGTH RANGE ============
    
    wavelengths = np.linspace(200, 750, 1000)
    absorption = np.zeros_like(wavelengths)


    # ============ BUILD THE ABSORPTION SPECTRA FROM PARSED DATA ============
    
    for peak_center, fosc, sigma in tddft_peaks:
        absorption += (fosc / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((wavelengths - peak_center) / sigma)**2)


    # ============ NORMALISE ABSORPTION AND CONVERT TO TRANSMITTANCE ============
    
    absorption_normalized = (absorption / np.max(absorption)) if np.max(absorption) > 0 else np.zeros_like(absorption)
    transmittance = 10 ** (-absorption_normalized)



    # ============ FIND PEAKS ============

    major_absorbance_peaks = find_major_peaks(absorption_normalized, wavelengths, min_prominence=0.01, min_height=0.01)
    major_transmittance_peaks = find_major_peaks(transmittance, wavelengths, min_prominence=0.01, min_height=0.01)    
    global_max_idx = np.argmax(transmittance)
    global_max_wl, global_max_trans = wavelengths[global_max_idx], transmittance[global_max_idx] 
    
    print("\nTop 5 Major Absorbance Peaks (λ_max values):\n" + "-" * 60)
    for i, (wl, abs_value) in enumerate(major_absorbance_peaks[:5]):
        print(f"λ_max {i+1}: {wl:.1f} nm (Abs = {abs_value:.3f})")

    print("\n\nTop 5 Major Transmittance Peaks:\n" + "-" * 60)
    for i, (wl, trans_value) in enumerate(major_transmittance_peaks[:5]):
        print(f"T_max {i+1}: {wl:.1f} nm (T = {trans_value:.3f} → {trans_value*100:.1f}%)")
    
    
    # ============ PREPARE STICK DATA ============
    
    stick_wavelengths = [p[0] for p in tddft_peaks]
    stick_fosc = [p[1] for p in tddft_peaks]
    max_fosc = max(stick_fosc) if stick_fosc else 1.0
    
    
    # ============ CREATE FIGURE ============
    
    nrows = 2 if show_transmittance else 1
    fig, axes = subplots(nrows, 1, figsize=(16, 12 if show_transmittance else 8))
    axes = [axes] if nrows == 1 else axes
    
    
    # ============ PLOT ABSORBANCE ============
    
    ax = axes[0]
    ax.plot(wavelengths, absorption_normalized, 'r-', label='Normalized Absorption', linewidth=3)
    
    
    # ============ ADD STICKS ============
    
    for wl, fosc in zip(stick_wavelengths, stick_fosc):
        if 200 <= wl <= 750:
            stick_height = (fosc / max_fosc) * stick_max_relative_height * max(absorption_normalized)
            ax.vlines(wl, 0, stick_height, colors=stick_colour, alpha=stick_alpha)
    
    add_spectrum_annotations(ax, major_absorbance_peaks, ['red', 'darkred', 'maroon', 'firebrick', 'crimson'], stick_wavelengths=stick_wavelengths, stick_fosc=stick_fosc, max_fosc=max_fosc, label_prefix='λ =')
    
    ax.set_ylabel('Normalized Absorbance (A)', size=18)
    ax.set_title(f'Simulated UV-Vis Spectrum from {path.basename(filename)} (σ = {sigma_value})')
    ax.grid(alpha=0.3)
    ax.axvline(380, color='gray', linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_ylim(-0.05, 1.05)
        
    # ============ PLOT TRANSMITTANCE IF REQUESTED ============

    if show_transmittance:
        ax = axes[1]
        ax.plot(wavelengths, transmittance, 'b-', label='Transmittance (T)', linewidth=3)
        
        add_spectrum_annotations(ax, major_transmittance_peaks, ['green', 'darkgreen', 'forestgreen', 'limegreen', 'seagreen'], stick_wavelengths=stick_wavelengths, stick_fosc=stick_fosc, max_fosc=max_fosc, label_prefix='λ =')
        
        ax.plot(global_max_wl, global_max_trans, 'mo', markersize=12, markeredgecolor='black', markeredgewidth=1, label='Highest Transmittance')
        ax.annotate(f'Max T = {global_max_wl:.1f} nm', xy=(global_max_wl, global_max_trans), xytext=(20, 20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.3', facecolor='pink', alpha=0.9, edgecolor='black'), fontsize=14)
        
        ax.set_xlabel('Wavelength (nm)', size=18)
        ax.set_ylabel('Transmittance (T)', size=18)
        ax.legend(loc='upper right', prop={'size': 18})
        ax.grid(alpha=0.3)
        ax.axvline(380, color='gray', linestyle='--', alpha=0.7)
        ax.tick_params(axis='both', labelsize=16)
        ax.set_ylim(-0.05, 1.05)
    else:
        axes[0].set_xlabel('Wavelength (nm)', size=18)
    
    tight_layout()
    
    
    # ============ SAVE FIGURE ============
    
    plot_type = "spectrum" if show_transmittance else "absorbance"
    output_filename = f"{path.splitext(filename)[0]}_{plot_type}_sigma{sigma_value}.png"
    savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n\nGraph saved as: {output_filename}")
    show()
    
    results = {'absorbance_peaks': major_absorbance_peaks, 'transmittance_peaks': major_transmittance_peaks, 'global_max_trans': (global_max_wl, global_max_trans), 'output_filename': output_filename, 'sigma_value': sigma_value}
    
    return results



def parse_tddft_file(filename, sigma_value):
    """ Parse TD-DFT data from a file.
    Returns list of [wavelength, fosc, sigma] for each transition. """
    
    tddft_peaks = []
    
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except Exception as e:
        print(f"Error reading file: {e}")
        return []
    
    data_started = False

    for line in lines:
        if "ABSORPTION SPECTRUM" in line or "-----" in line or not line.strip():
            continue
        
        if "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:
            break
            
        if data_started or match(r'\s*\d+-\d+[A-Z]\s+->\s+\d+-\d+[A-Z]', line):
            data_started = True
            
            m = match(r'\s*[\d\-A-Z]+\s+->\s+[\d\-A-Z]+\s+([\d\.]+)\s+[\d\.]+\s+([\d\.]+)\s+([\d\.]+)', line)
            if m:
                wavelength = float(m.group(2))
                fosc = float(m.group(3))
                                
                if 200 <= wavelength <= 780:
                    tddft_peaks.append([wavelength, fosc, sigma_value])
    
    return tddft_peaks



def wavelength_to_rgb(wavelength, gamma=0.8):
    """ Convert wavelength to RGB colour.
    Return RGB values as floats between 0 and 1: (R, G, B). """
    
    wavelength = max(380, min(750, float(wavelength)))
    
    if wavelength < 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / 60
        R, G, B = ((-(wavelength - 440) / 60) * attenuation) ** gamma, 0.0, attenuation ** gamma
    elif wavelength < 490:
        R, G, B = 0.0, ((wavelength - 440) / (490 - 440)) ** gamma, 1.0
    elif wavelength < 510:
        R, G, B = 0.0, 1.0, (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength < 580:
        R, G, B = ((wavelength - 510) / (580 - 510)) ** gamma, 1.0, 0.0
    elif wavelength < 645:
        R, G, B = 1.0, (-(wavelength - 645) / (645 - 580)) ** gamma, 0.0
    elif wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R, G, B = (1.0 * attenuation) ** gamma, 0.0, 0.0
    else:
        R, G, B = 0.0, 0.0, 0.0
        
    return (R, G, B)



# Example usage:
if __name__ == "__main__":
    filename = "Org1b+2a.out"                       # Change this to your file
    results = generate_spectrum(filename)