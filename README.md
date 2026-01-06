<<<<<<< HEAD
# React + Vite

This template provides a minimal setup to get React working in Vite with HMR and some ESLint rules.

Currently, two official plugins are available:

- [@vitejs/plugin-react](https://github.com/vitejs/vite-plugin-react/blob/main/packages/plugin-react) uses [Babel](https://babeljs.io/) (or [oxc](https://oxc.rs) when used in [rolldown-vite](https://vite.dev/guide/rolldown)) for Fast Refresh
- [@vitejs/plugin-react-swc](https://github.com/vitejs/vite-plugin-react/blob/main/packages/plugin-react-swc) uses [SWC](https://swc.rs/) for Fast Refresh

## React Compiler

The React Compiler is not enabled on this template because of its impact on dev & build performances. To add it, see [this documentation](https://react.dev/learn/react-compiler/installation).

## Expanding the ESLint configuration

If you are developing a production application, we recommend using TypeScript with type-aware lint rules enabled. Check out the [TS template](https://github.com/vitejs/vite/tree/main/packages/create-vite/template-react-ts) for information on how to integrate TypeScript and [`typescript-eslint`](https://typescript-eslint.io) in your project.
=======
# ACOUSTIC WORKBENCH  
### 9-Parameter In-Ear Acoustic Solver (v4.3)

ACOUSTIC WORKBENCH is an interactive, browser-based acoustic simulation and analysis tool for **in-ear headphone design**, built on a **physics-based Equivalent Circuit Model (ECM)**.  
The system enables real-time exploration of cavity, leakage, and damping parameters, and provides **frequency response, impulse response, and perceptual timbre analysis** within a unified UI.

This project is designed for **acoustics research, headphone prototyping, and educational use**, and serves as the frontend foundation for future **PINN-based inverse design and optimization workflows**.

---

## ✨ Key Features

### 🔧 9-Parameter Physical Model (ECM)
The solver is based on a lumped-parameter electro-acoustic equivalent circuit, with **9 tunable physical parameters**:

| Category | Parameter |
|--------|-----------|
| Rear Cavity | Rear Volume (Var2) |
| Leakage Path | Leak Length (a), Leak Width (b) |
| Tube / Nozzle | Tube Length (ft1), Tube Radius (fat1) |
| Leakage Dynamics | Leak Mass (mleakage), Leak Resistance (Rleakage) |
| Damping | Damping Resistance (Rp1), Damping Mass (mp1) |

All parameters are adjusted via real-time sliders and mapped directly to physical units.

---

### 📈 Frequency Response (SPL)
- Log-scale SPL from **20 Hz – 20 kHz**
- Real-time ECM simulation
- **Harman Target** curve with automatic level alignment (1 kHz reference)
- Optional **measured WAV import** for comparison

---

### ⏱ Impulse Response (IR)
- Time-domain impulse response reconstructed from simulated magnitude
- Normalized peak alignment
- Direct comparison with measured IR (if loaded)

---

### 🎧 Timbre Analysis (Perceptual Layer)
The system performs band-wise deviation analysis relative to the Harman target:

- Sub-Bass
- Bass
- Mid
- Treble
- Air

Outputs include:
- **Qualitative timbre tags** (e.g. *Bright*, *Lack Sub Bass*, *Lack of Air*)
- **Quantitative spectrum deviation bars** for fast tuning feedback

---

### ⚙️ Auto Fit (Heuristic Optimization)
- Random-target generation for testing inverse behavior
- Automatic parameter fitting via stochastic hill-climbing
- Designed for future replacement with **Physics-Informed Neural Networks (PINN)**

---

### 🚀 Performance Optimizations
- **150 ms debounced parameter updates** to prevent excessive recomputation
- Immediate (non-debounced) updates during Auto Fit
- Client-side only, no backend required

---

## 🖥 UI Overview

- Dark, research-oriented UI
- Left panel: physical parameters
- Main panel:
  - SPL response
  - Impulse response
  - Timbre analysis & spectrum deviation
- View switching between **Home / Circuit**

> The UI is designed to resemble professional acoustic and measurement software rather than consumer EQ tools.

---

## 🛠 Tech Stack

- **Frontend Framework**: React (Vite)
- **Visualization**: Recharts
- **Icons**: lucide-react
- **Styling**: Tailwind CSS
- **Math Core**: Custom complex arithmetic, ECM solver, FFT-based IR reconstruction

---

## 📂 Project Structure

```text
frontend/
├─ src/
│  ├─ InEarDesigner.jsx     # Main application logic & UI
│  ├─ App.jsx
│  ├─ main.jsx
│  └─ index.css             # Tailwind + custom styles
├─ public/
├─ index.html
├─ vite.config.js
├─ tailwind.config.js
├─ package.json
└─ README.md
>>>>>>> a84d9c9 (docs: add project README)
