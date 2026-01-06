import React, { useState, useEffect, useRef, useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Label, ReferenceLine, Legend, Scatter, ComposedChart } from 'recharts';
import { Activity, Settings, Save, FileText, BarChart2, Sliders, Zap, Waves, CircuitBoard, RotateCcw, Calculator, ChevronDown, ChevronUp, Upload, Play, Shuffle } from 'lucide-react';

// --- 1. Math Kernel (Complex & FFT) ---

class Complex {
  constructor(r, i) { this.r = r; this.i = i; }
  static add(a, b) { return new Complex(a.r + b.r, a.i + b.i); }
  static sub(a, b) { return new Complex(a.r - b.r, a.i - b.i); }
  static mul(a, b) { return new Complex(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
  static div(a, b) {
    const den = b.r * b.r + b.i * b.i;
    if (den === 0) return new Complex(0, 0);
    return new Complex((a.r * b.r + a.i * b.i) / den, (a.i * b.r - a.r * b.i) / den);
  }
  static abs(a) { return Math.sqrt(a.r * a.r + a.i * a.i); }
  static fromReal(r) { return new Complex(r, 0); }
  static fromImag(i) { return new Complex(0, i); }
}

const solveComplexMatrix = (A, B) => {
  const n = B.length;
  const M = A.map(row => [...row]);
  const x = [...B]; 

  for (let i = 0; i < n; i++) {
    let maxRow = i;
    let maxVal = Complex.abs(M[i][i]);
    for (let k = i + 1; k < n; k++) {
      if (Complex.abs(M[k][i]) > maxVal) {
        maxVal = Complex.abs(M[k][i]);
        maxRow = k;
      }
    }
    [M[i], M[maxRow]] = [M[maxRow], M[i]];
    [x[i], x[maxRow]] = [x[maxRow], x[i]];

    for (let k = i + 1; k < n; k++) {
      const factor = Complex.div(M[k][i], M[i][i]);
      x[k] = Complex.sub(x[k], Complex.mul(factor, x[i]));
      for (let j = i; j < n; j++) {
        M[k][j] = Complex.sub(M[k][j], Complex.mul(factor, M[i][j]));
      }
    }
  }

  const res = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let sum = new Complex(0, 0);
    for (let j = i + 1; j < n; j++) {
      sum = Complex.add(sum, Complex.mul(M[i][j], res[j]));
    }
    res[i] = Complex.div(Complex.sub(x[i], sum), M[i][i]);
  }
  return res;
};

// FFT
const fft = (data) => {
  const N = data.length;
  if (N <= 1) return data;
  const half = N / 2;
  const even = new Array(half);
  const odd = new Array(half);
  for (let i = 0; i < half; i++) {
    even[i] = data[2 * i];
    odd[i] = data[2 * i + 1];
  }
  const evenResult = fft(even);
  const oddResult = fft(odd);
  const result = new Array(N);
  for (let k = 0; k < half; k++) {
    const angle = -2 * Math.PI * k / N;
    const w = new Complex(Math.cos(angle), Math.sin(angle));
    const t = Complex.mul(w, oddResult[k]);
    result[k] = Complex.add(evenResult[k], t);
    result[k + half] = Complex.sub(evenResult[k], t);
  }
  return result;
};

// --- 2. Physics Constants ---
const PI = Math.PI;
const CONST = {
  rho0: 1.21, c: 343.0,
  Re: 15.818, Le: 3.5e-5,
  Rm: 0.75619, mm: 2.0708e-4, Cm: 1.3502e-05, Bl: 1.6697,
  eg: 64.0, AD: 0.5e-3, mu: 18.6e-6,
  Bi: 3 * PI / 8, Br: 0.5,
  Vaf: 46.13e-8, 
  Var: 20e-8,    
  at: 2.5e-3, lt: 17e-3,
  ma1: 82.9, Ca1: 0.943e-12, Ra2: 50.6e6, ma2: 9.4e3, Ca2: 1.9e-12,
  ma3: 130.3, Ca3: 1.492e-12, Ra4: 31.1e6, ma4: 983.8, Ca4: 2.1e-12,
  ma5: 133.4, Ca5: 1.517e-12,
  at2: 2.8e-3, lt2: 2.8e-3,
  af1: (2.5e-3)/2, af2: (1.57e-3)/2, lf2: 2e-3, ar1: (1.5e-3)/2
};

const HARMAN_TARGET = [
  { f: 20, spl: 125 }, { f: 30, spl: 126 }, { f: 40, spl: 126 }, { f: 60, spl: 125 },
  { f: 100, spl: 123 }, { f: 200, spl: 119 }, { f: 300, spl: 117 }, { f: 400, spl: 116 },
  { f: 1000, spl: 118 }, { f: 2000, spl: 124 }, { f: 3000, spl: 128 }, { f: 4000, spl: 126 },
  { f: 6000, spl: 120 }, { f: 8000, spl: 115 }, { f: 10000, spl: 110 }, { f: 16000, spl: 95 }, { f: 20000, spl: 85 }
];

const getTargetSPL = (freq) => {
  for (let i = 0; i < HARMAN_TARGET.length - 1; i++) {
    if (freq >= HARMAN_TARGET[i].f && freq <= HARMAN_TARGET[i+1].f) {
      const range = HARMAN_TARGET[i+1].f - HARMAN_TARGET[i].f;
      const dist = freq - HARMAN_TARGET[i].f;
      const interpFactor = dist / range; 
      return HARMAN_TARGET[i].spl + interpFactor * (HARMAN_TARGET[i+1].spl - HARMAN_TARGET[i].spl);
    }
  }
  return 100;
};

// Fixed terms
const AF1 = PI * CONST.af1**2;
const RAF1 = 8 * CONST.mu * CONST.lt / (PI * CONST.af1**4); 
const MAF1 = 4 * CONST.rho0 * CONST.af1 / (4 * AF1); 
const AF2 = PI * CONST.af2**2;
const RAF2 = 8 * CONST.mu * CONST.lt / (PI * CONST.af2**4);
const MAF2 = 4 * CONST.rho0 * CONST.af2 / (4 * AF2);
const RAF12 = RAF1 + RAF2;
const MAF12 = MAF1 + MAF2;
const AT = PI * CONST.at**2;
const MAT = CONST.rho0 * CONST.lt / (2 * AT);
const CAT = AT * CONST.lt / (CONST.rho0 * CONST.c**2);
const CAF = CONST.Vaf / (CONST.rho0 * CONST.c**2);
const AR1 = PI * CONST.ar1**2;
const RAR1 = 8 * CONST.mu * 1.5 / (PI * CONST.ar1**4);
const MAR1 = 4 * CONST.rho0 * CONST.ar1 / (4 * AR1);
const AT2 = PI * CONST.at2**2;
const MFT2 = CONST.rho0 * CONST.lt2 / (2 * AT2);
const CFT2 = (AT2 * CONST.lt2) / (CONST.rho0 * CONST.c**2);

// --- 3. Engines ---

const calculateECM = (params) => {
  // 9 Parameters
  const { Var2, a, b, ft1, fat1, mleakage, Rleakage, Rp1, mp1 } = params;

  const Car2 = Var2 / (CONST.rho0 * CONST.c**2);
  const Af = a * b; 
  const P = 2 * (a + b); 
  const rh = Af / (P + 1e-20);
  const Rar_2 = 12.0 * CONST.mu * CONST.lf2 / (a * b**3 + 1e-20);
  const Mar_2 = CONST.rho0 * rh / (Af + 1e-20);
  const Rarad = (CONST.Bi**2) * CONST.rho0 * CONST.c / (CONST.Br * 4.0 * Af + 1e-20);
  const Marad = CONST.Bi * CONST.rho0 * rh / (4.0 * Af + 1e-20);

  const fAt1 = PI * fat1**2;
  const mft1 = CONST.rho0 * ft1 / (2.0 * fAt1 + 1e-20);
  const Cft1 = fAt1 * ft1 / (CONST.rho0 * CONST.c**2 + 1e-20);

  const rawCurve = [];
  const steps = 100; const fMin = 20; const fMax = 20000;
  const one = new Complex(1, 0);
  
  for (let i = 0; i <= steps; i++) {
    const f = fMin * Math.pow(fMax / fMin, i / steps);
    const w = 2 * PI * f; const jw = new Complex(0, w);

    // Impedance calculation
    const Z_electrical = Complex.add(Complex.fromReal(CONST.Re), Complex.mul(jw, Complex.fromReal(CONST.Le)));
    const Zeb = Complex.div(Complex.fromReal(CONST.Bl**2), Complex.mul(Complex.fromReal(CONST.AD**2), Z_electrical));
    const Z_mech_series = Complex.add(Complex.add(Complex.fromReal(CONST.Rm), Complex.mul(jw, Complex.fromReal(CONST.mm))), Complex.div(one, Complex.mul(jw, Complex.fromReal(CONST.Cm))));
    const Zm = Complex.div(Z_mech_series, Complex.fromReal(CONST.AD**2));
    const Z1 = Complex.add(Zeb, Zm); 

    const CAR_val = CONST.Var / (CONST.rho0 * CONST.c**2);
    const Z2 = Complex.div(one, Complex.mul(jw, Complex.fromReal(CAR_val)));
    const Z3 = Complex.add(Complex.fromReal(RAR1), Complex.mul(jw, Complex.fromReal(MAR1)));
    const Z4 = Complex.div(one, Complex.mul(jw, Complex.fromReal(Car2)));

    const jwMarad = Complex.mul(jw, Complex.fromReal(Marad));
    const Z_rad = Complex.div(Complex.mul(Complex.fromReal(Rarad), jwMarad), Complex.add(Complex.fromReal(Rarad), jwMarad));
    const Z5 = Complex.add(Complex.add(Complex.fromReal(Rar_2), Complex.mul(jw, Complex.fromReal(Mar_2))), Z_rad);

    const Z6 = Complex.div(one, Complex.mul(jw, Complex.fromReal(CAF)));
    const Z7 = Complex.add(Complex.add(Complex.fromReal(RAF12), Complex.mul(jw, Complex.fromReal(MAF12))), Complex.mul(jw, Complex.fromReal(mft1 / 2.0)));
    const Z8 = Complex.div(one, Complex.mul(jw, Complex.fromReal(Cft1)));
    const Z9 = Complex.mul(jw, Complex.fromReal(mft1 / 2.0));
    const Z10 = Complex.add(Complex.fromReal(Rleakage), Complex.mul(jw, Complex.fromReal(mleakage)));

    const Z11 = Complex.mul(jw, Complex.fromReal(MAT / 2.0));
    const Z12 = Complex.div(one, Complex.mul(jw, Complex.fromReal(CAT)));
    const Z13 = Complex.add(Complex.mul(jw, Complex.fromReal(MAT / 2.0)), Complex.mul(jw, Complex.fromReal(CONST.ma1)));
    const Z14 = Complex.div(one, Complex.mul(jw, Complex.fromReal(CONST.Ca1)));
    const Z15 = Complex.add(Complex.add(Complex.fromReal(CONST.Ra2), Complex.mul(jw, Complex.fromReal(CONST.ma2))), Complex.div(one, Complex.mul(jw, Complex.fromReal(CONST.Ca2))));
    const Z16 = Complex.mul(jw, Complex.fromReal(CONST.ma3));
    const Z17 = Complex.div(one, Complex.mul(jw, Complex.fromReal(CONST.Ca3)));
    const Z18 = Complex.add(Complex.add(Complex.fromReal(CONST.Ra4), Complex.mul(jw, Complex.fromReal(CONST.ma4))), Complex.div(one, Complex.mul(jw, Complex.fromReal(CONST.Ca4))));
    const Z19 = Complex.mul(jw, Complex.fromReal(CONST.ma5));
    const Z20 = Complex.div(one, Complex.mul(jw, Complex.fromReal(CONST.Ca5)));

    const Z21 = Complex.mul(jw, Complex.fromReal(MFT2 / 2.0));
    const Z22 = Complex.div(one, Complex.mul(jw, Complex.fromReal(CFT2)));
    const Z23 = Complex.mul(jw, Complex.fromReal(MFT2 / 2.0));
    
    const Zmp1 = Complex.mul(jw, Complex.fromReal(mp1));
    const Rp1_c = Complex.fromReal(Rp1);
    const ZP = Complex.div(Complex.mul(Rp1_c, Zmp1), Complex.add(Rp1_c, Zmp1));

    const zero = new Complex(0, 0);
    const A = Array(12).fill(0).map(() => Array(12).fill(zero));
    
    A[0][0] = Complex.add(Complex.add(Z1, Z2), Z6); A[0][1] = Complex.mul(Z2, new Complex(-1, 0)); A[0][3] = Complex.mul(Z6, new Complex(-1, 0));
    A[1][0] = Complex.mul(Z2, new Complex(-1, 0)); A[1][1] = Complex.add(Complex.add(Z2, Z3), Z4); A[1][2] = Complex.mul(Z4, new Complex(-1, 0));
    A[2][1] = Complex.mul(Z4, new Complex(-1, 0)); A[2][2] = Complex.add(Z4, Z5);
    A[3][0] = Complex.mul(Z6, new Complex(-1, 0)); A[3][3] = Complex.add(Complex.add(Z6, Z7), Z8); A[3][4] = Complex.mul(Z8, new Complex(-1, 0));
    A[4][3] = Complex.mul(Z8, new Complex(-1, 0)); A[4][4] = Complex.add(Complex.add(Z8, Z9), Z10); A[4][5] = Complex.mul(Z10, new Complex(-1, 0));
    A[5][4] = Complex.mul(Z10, new Complex(-1, 0)); A[5][5] = Complex.add(Complex.add(Z10, Z21), Z22); A[5][6] = Complex.mul(Z22, new Complex(-1, 0));
    A[6][5] = Complex.mul(Z22, new Complex(-1, 0)); A[6][6] = Complex.add(Complex.add(Complex.add(Z22, Z23), ZP), Complex.add(Z11, Z12)); A[6][7] = Complex.mul(Z12, new Complex(-1, 0));
    A[7][6] = Complex.mul(Z12, new Complex(-1, 0)); A[7][7] = Complex.add(Complex.add(Z12, Z13), Z14); A[7][8] = Complex.mul(Z14, new Complex(-1, 0));
    A[8][7] = Complex.mul(Z14, new Complex(-1, 0)); A[8][8] = Complex.add(Z14, Z15); A[8][9] = Complex.mul(Z15, new Complex(-1, 0));
    A[9][8] = Complex.mul(Z15, new Complex(-1, 0)); A[9][9] = Complex.add(Complex.add(Z15, Z16), Z17); A[9][10] = Complex.mul(Z17, new Complex(-1, 0));
    A[10][9] = Complex.mul(Z17, new Complex(-1, 0)); A[10][10] = Complex.add(Z17, Z18); A[10][11] = Complex.mul(Z18, new Complex(-1, 0));
    A[11][10] = Complex.mul(Z18, new Complex(-1, 0)); A[11][11] = Complex.add(Complex.add(Z18, Z19), Z20);

    const Vin_num = Complex.fromReal(CONST.eg * CONST.Bl / CONST.AD);
    const Vin = Complex.div(Vin_num, Z_electrical);
    const B_vec = Array(12).fill(zero);
    B_vec[0] = Vin;

    const I = solveComplexMatrix(A, B_vec);
    const H = Complex.mul(I[11], Z20);
    const mag = Complex.abs(H);
    const spl = 20 * Math.log10(mag + 1e-20);
    
    // Store HARMAN Target as "targetRaw" for later alignment
    const targetSpl = getTargetSPL(f);

    rawCurve.push({
      frequency: Math.round(f),
      spl: parseFloat(spl.toFixed(2)),
      targetRaw: targetSpl,
      mag: mag 
    });
  }

  const tubeResonanceFreq = CONST.c / (4 * (ft1 + 0.6 * 2 * fat1 + 0.015)); 
  const cutoffFreq = 1 / Math.log10(Rleakage + 1) * 100; 

  return { curveData: rawCurve, metrics: { tubeResonanceFreq, cutoffFreq } }; 
};

const calculateIR = (curveData) => {
  const N = 1024; const fs = 48000;
  const mags = new Array(N).fill(0);
  
  // Interpolate to linear grid for FFT
  for (let k = 0; k < N/2; k++) {
    const f = k * fs / N;
    if (f < 20) continue; 
    const pt = curveData.reduce((prev, curr) => Math.abs(curr.frequency - f) < Math.abs(prev.frequency - f) ? curr : prev);
    mags[k] = pt.mag;
  }

  const rawIR = new Array(N).fill(0);
  for (let n = 0; n < N; n++) { 
    let sum = 0;
    for (let k = 1; k < N/2; k++) {
      const w = 2 * PI * k * n / N; sum += mags[k] * Math.cos(w);
    }
    rawIR[n] = sum;
  }
  
  let maxVal = Math.abs(rawIR[0]); // Peak at 0 for zero phase
  const viewDurationMs = 3.0; 
  const samples = Math.floor(viewDurationMs / 1000 * fs);
  const irData = [];

  for (let i = 0; i < samples; i++) {
    const t = (i / fs * 1000); 
    irData.push({ time: parseFloat(t.toFixed(3)), amplitude: rawIR[i] / (maxVal || 1) });
  }
  
  return { irData, peakPoint: { time: 0, amplitude: 1.0 } };
};

// --- Analysis Logic ---
const analyzeTimbre = (data, thresholds) => {
  const bands = {
    subBass: { min: 20, max: 60, avg: 0, count: 0 },
    bass: { min: 60, max: 250, avg: 0, count: 0 },
    lowMid: { min: 250, max: 500, avg: 0, count: 0 },
    mid: { min: 500, max: 2000, avg: 0, count: 0 },
    upperMid: { min: 2000, max: 5000, avg: 0, count: 0 },
    treble: { min: 5000, max: 10000, avg: 0, count: 0 },
    air: { min: 10000, max: 20000, avg: 0, count: 0 }
  };

  data.forEach(p => {
    // Check if 'delta' exists (it might not if target isn't aligned yet)
    // We compute delta on the fly if needed
    const delta = (p.delta !== undefined) ? p.delta : (p.spl - p.target);
    
    for (const key in bands) {
      if (p.frequency >= bands[key].min && p.frequency < bands[key].max) {
        bands[key].avg += delta;
        bands[key].count++;
      }
    }
  });

  const tags = [];
  const analysis = {};
  const { sensitivity, sharpness } = thresholds;

  for (const key in bands) {
    const count = bands[key].count || 1;
    const avgDelta = bands[key].avg / count;
    analysis[key] = avgDelta;

    if (avgDelta > sensitivity) {
      if (key === 'subBass') tags.push("Strong Sub Bass");
      if (key === 'bass') tags.push("Bass Heavy");
      if (key === 'lowMid' || key === 'mid') tags.push("Mid Forward");
      if (key === 'upperMid' || key === 'treble') tags.push("Bright");
      if (key === 'air') tags.push("Airy");
    } else if (avgDelta < -sensitivity) {
      if (key === 'subBass') tags.push("Lack Sub Bass");
      if (key === 'bass') tags.push("Bass Light");
      if (key === 'lowMid' || key === 'mid') tags.push("Recessed Mid");
      if (key === 'upperMid') tags.push("Smooth");
      if (key === 'treble') tags.push("Dark Treble");
      if (key === 'air') tags.push("Lack of Air");
    }
  }

  if ((analysis.upperMid || 0) > sharpness) tags.push("Harsh/Sharp");
  if ((analysis.bass || 0) > (sensitivity * 1.5)) tags.push("Boomy");

  return { tags: Array.from(new Set(tags)), analysis };
};

// --- Helper: Format & Scaling ---
const Scale = {
  // Mapping UI sliders to Physical units
  ranges: {
    Var2:     { min: 1.48e-4, max: 6.00e-4, step: 0.01e-4, unit: 'm³', label: 'Rear Vol' },
    a:        { min: 0.70e-3, max: 1.50e-3, step: 0.05e-3, unit: 'm',  label: 'Leak L' },
    b:        { min: 0.20e-3, max: 0.80e-3, step: 0.05e-3, unit: 'm',  label: 'Leak W' },
    ft1:      { min: 8.00e-3, max: 15.0e-3, step: 0.10e-3, unit: 'm',  label: 'Tube L' },
    fat1:     { min: 1.25e-3, max: 3.00e-3, step: 0.05e-3, unit: 'm',  label: 'Tube R' },
    mleakage: { min: 1.00e4,  max: 2.00e5,  step: 0.10e4,  unit: 'kg/m4', label: 'Leak Mass' },
    Rleakage: { min: 0.05e9,  max: 1.50e10, step: 0.05e9,  unit: 'Ω',  label: 'Leak Res' },
    Rp1:      { min: 2.50e6,  max: 10.0e6,  step: 0.10e6,  unit: 'Ω',  label: 'Damp R' },
    mp1:      { min: 40.0,    max: 160.0,   step: 1.0,     unit: '',   label: 'Damp M' },
  }
};

const formatNumber = (num) => {
  if (Math.abs(num) < 1e-6 || Math.abs(num) > 1e4) return num.toExponential(2);
  return num.toPrecision(3);
}

// --- Components ---

const SliderControl = ({ name, value, config, onChange }) => (
  <div className="mb-4 group">
    <div className="flex justify-between items-center mb-1">
      <label className="text-xs text-slate-400 font-medium group-hover:text-cyan-400 transition-colors flex items-center gap-1">
        {config.label} <span className="text-[10px] opacity-50">({name})</span>
      </label>
      <span className="text-xs font-mono text-cyan-400 bg-cyan-950/30 px-1.5 py-0.5 rounded border border-cyan-900/50">
        {formatNumber(value)} <span className="text-slate-500">{config.unit}</span>
      </span>
    </div>
    <div className="flex items-center gap-2">
      <input 
        type="range" min={config.min} max={config.max} step={config.step} value={value} 
        onChange={(e) => onChange(name, parseFloat(e.target.value))}
        className="w-full h-1 bg-slate-700 rounded-lg appearance-none cursor-pointer accent-cyan-500 hover:accent-cyan-400"
      />
    </div>
  </div>
);

const CircuitDiagram = () => (
  <div className="w-full h-full flex items-center justify-center bg-slate-900 p-8 select-none text-slate-500">
    <div className="text-center">
      <CircuitBoard size={64} className="mx-auto mb-4 opacity-50"/>
      <p>ECM Circuit Topology (Z1...Z23)</p>
      <p className="text-xs opacity-60">Visual representation same as v3.8</p>
    </div>
  </div>
);

export default function InEarDesigner() {
  const [params, setParams] = useState({
    Var2: 2.0e-4, a: 1.0e-3, b: 0.5e-3, ft1: 11e-3, fat1: 2.0e-3,
    mleakage: 5e4, Rleakage: 1e9, Rp1: 5e6, mp1: 80.0
  });

  const [measurement, setMeasurement] = useState(null); // { spl: [], ir: [] }
  const [analysisThresholds, setAnalysisThresholds] = useState({ sensitivity: 2.5, sharpness: 5.0 });
  const [viewMode, setViewMode] = useState('spl'); 
  const [chartData, setChartData] = useState([]);
  const [irData, setIrData] = useState([]);
  const [irPeak, setIrPeak] = useState(null);
  const [timbre, setTimbre] = useState({ tags: [], analysis: {} });
  const [isOptimizing, setIsOptimizing] = useState(false);
  const fileInputRef = useRef(null);
  const optimRef = useRef(null);

  // --- Core Loop ---
  useEffect(() => {
    // 1. Calculate Simulation
    const { curveData: rawSimData } = calculateECM(params);
    const { irData: ir, peakPoint } = calculateIR(rawSimData);

    // 2. Prepare Display Data & Align Target
    const refFreq = 1000;
    const simPt = rawSimData.reduce((p, c) => Math.abs(c.frequency - refFreq) < Math.abs(p.frequency - refFreq) ? c : p);
    const harmanOffset = simPt.spl - simPt.targetRaw;

    let displayData = rawSimData.map(p => {
      const alignedTarget = p.targetRaw + harmanOffset;
      let measSpl = null;
      if (measurement) {
        const mPt = measurement.spl.find(mp => Math.abs(mp.frequency - p.frequency) < 5);
        if (mPt) measSpl = mPt.spl;
      }
      return {
        ...p,
        target: parseFloat(alignedTarget.toFixed(2)), 
        meas: measSpl, 
        delta: parseFloat((p.spl - alignedTarget).toFixed(2)) 
      };
    });

    // 3. Prepare IR Data with Measurement Comparison
    let displayIrData = ir.map(p => ({...p, measAmp: null}));
    if (measurement && measurement.ir) {
      // Merge Sim IR and Meas IR
      displayIrData = ir.map((p, i) => ({
        ...p,
        measAmp: measurement.ir[i] ? measurement.ir[i].amplitude : null
      }));
    }

    setChartData(displayData);
    setIrData(displayIrData);
    setIrPeak(peakPoint);
    setTimbre(analyzeTimbre(displayData, analysisThresholds));
  }, [params, measurement, analysisThresholds]);

  // --- Handlers ---

  const handleParamChange = (name, val) => {
    setParams(prev => ({...prev, [name]: val}));
  };

  const handleFileUpload = async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    try {
      const arrayBuffer = await file.arrayBuffer();
      const audioCtx = new (window.AudioContext || window.webkitAudioContext)();
      const audioBuffer = await audioCtx.decodeAudioData(arrayBuffer);
      const channelData = audioBuffer.getChannelData(0); 
      
      const fftSize = 2048;
      const input = new Array(fftSize).fill(0);
      for(let i=0; i<Math.min(channelData.length, fftSize); i++) input[i] = new Complex(channelData[i], 0);
      const output = fft(input);
      
      const measSpl = [];
      const fs = audioBuffer.sampleRate;
      for(let i=0; i<fftSize/2; i++) {
        const freq = i * fs / fftSize;
        if(freq < 20 || freq > 20000) continue;
        const mag = Complex.abs(output[i]);
        const db = 20 * Math.log10(mag + 1e-10); 
        measSpl.push({ frequency: freq, spl: db });
      }
      
      // Also process IR for display (Raw time data)
      const maxVal = Math.max(...Array.from(channelData).map(Math.abs));
      const samples = Math.floor(3.0 / 1000 * fs); // 3ms
      let peakIdx = 0; let peakV = 0;
      for(let i=0; i<channelData.length; i++) { if(Math.abs(channelData[i])>peakV){ peakV=Math.abs(channelData[i]); peakIdx=i;}}
      
      const measIr = [];
      for(let i=0; i<samples; i++) {
        const idx = peakIdx + i; 
        if(idx < channelData.length) {
           measIr.push({ time: parseFloat((i/fs*1000).toFixed(3)), amplitude: channelData[idx]/maxVal });
        }
      }

      setMeasurement({ spl: measSpl, ir: measIr });
    } catch (err) {
      console.error("WAV load failed", err);
      alert("Failed to load WAV.");
    }
  };

  const generateRandomTarget = () => {
    const randParams = {};
    Object.keys(Scale.ranges).forEach(k => {
      const r = Scale.ranges[k];
      // Increase variance significantly as requested
      // For Rleakage (Log scale behavior), standard rand is fine for linear sliders?
      // No, ranges are linear here.
      // ft1: Tube Length. min 8e-3, max 15e-3. 
      // Let's use full range uniformly.
      randParams[k] = r.min + Math.random() * (r.max - r.min);
    });
    
    // Calculate full curve and IR for target
    const { curveData } = calculateECM(randParams);
    const { irData } = calculateIR(curveData);
    
    setMeasurement({ 
      spl: curveData.map(p => ({ frequency: p.frequency, spl: p.spl })),
      ir: irData
    });
  };

  const runOptimization = () => {
    if (isOptimizing) {
      setIsOptimizing(false);
      cancelAnimationFrame(optimRef.current);
      return;
    }
    if (!measurement) return;

    setIsOptimizing(true);
    
    const step = () => {
      setParams(current => {
        const next = { ...current };
        const keys = Object.keys(Scale.ranges);
        const k = keys[Math.floor(Math.random() * keys.length)];
        const r = Scale.ranges[k];
        
        const range = r.max - r.min;
        const delta = range * 0.02 * (Math.random() > 0.5 ? 1 : -1); // Bigger step
        
        const trialParams = { ...current, [k]: current[k] + delta };
        if (trialParams[k] < r.min) trialParams[k] = r.min;
        if (trialParams[k] > r.max) trialParams[k] = r.max;

        const getLoss = (p) => {
           const res = calculateECM(p);
           const simMean = res.curveData.reduce((s,x)=>s+x.spl,0)/res.curveData.length;
           const measMean = measurement.spl.reduce((s,x)=>s+x.spl,0)/measurement.spl.length;
           const offset = measMean - simMean;

           let err = 0;
           let count = 0;
           res.curveData.forEach(pt => {
             const m = measurement.spl.find(mp => Math.abs(mp.frequency - pt.frequency) < 10);
             if(m) {
                err += (pt.spl - (m.spl - offset)) ** 2;
                count++;
             }
           });
           return count > 0 ? err/count : 1e9;
        };

        const l1 = getLoss(current);
        const l2 = getLoss(trialParams);

        if (l2 < l1) return trialParams;
        return current;
      });
      
      optimRef.current = requestAnimationFrame(step);
    };
    step();
  };

  useEffect(() => {
    if(!isOptimizing && optimRef.current) cancelAnimationFrame(optimRef.current);
  }, [isOptimizing]);


  return (
    <div className="flex flex-col h-screen bg-slate-950 text-slate-300 font-sans overflow-hidden">
      
      <header className="h-14 bg-slate-900 border-b border-slate-800 flex items-center px-6 justify-between shrink-0 shadow-sm z-10">
        <div className="flex items-center gap-3">
          <div className="w-8 h-8 bg-cyan-600 rounded flex items-center justify-center shadow-lg shadow-cyan-900/20">
            <Activity className="text-white w-5 h-5" />
          </div>
          <div>
            <h1 className="text-lg font-bold text-slate-100 tracking-tight">ACOUSTIC<span className="text-cyan-500">WORKBENCH</span></h1>
            <p className="text-[10px] text-slate-500 font-medium uppercase tracking-wider">9-Param Solver v4.3</p>
          </div>
        </div>
        <div className="flex items-center gap-4">
          <div className="flex gap-1 bg-slate-800 p-1 rounded-lg text-sm font-medium">
            <button onClick={() => setViewMode('spl')} className={`px-3 py-1.5 rounded-md flex items-center gap-2 transition-all ${viewMode === 'spl' ? 'bg-cyan-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200 hover:bg-slate-700'}`}><Waves size={16}/> Home</button>
            <button onClick={() => setViewMode('circuit')} className={`px-3 py-1.5 rounded-md flex items-center gap-2 transition-all ${viewMode === 'circuit' ? 'bg-cyan-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200 hover:bg-slate-700'}`}><CircuitBoard size={16}/> Circuit</button>
          </div>
        </div>
      </header>

      <div className="flex-1 flex overflow-hidden">
        {/* Left Sidebar: 9 Params */}
        <div className="w-80 bg-slate-900 border-r border-slate-800 flex flex-col shrink-0 overflow-y-auto custom-scrollbar">
          <div className="p-6 border-b border-slate-800">
            <div className="flex items-center justify-between mb-6">
              <div className="flex items-center gap-2 text-cyan-500">
                <Sliders className="w-5 h-5" />
                <h3 className="text-sm font-bold uppercase tracking-wider">Parameters (9)</h3>
              </div>
            </div>
            
            {Object.keys(Scale.ranges).map(key => (
              <SliderControl 
                key={key} 
                name={key} 
                value={params[key]} 
                config={Scale.ranges[key]} 
                onChange={handleParamChange} 
              />
            ))}
          </div>
        </div>

        {/* Main Content */}
        <div className="flex-1 flex flex-col min-w-0 bg-slate-950 relative">
          
          {/* Action Overlay (Only in Home View) */}
          {viewMode === 'spl' && (
            <div className="absolute top-4 right-6 z-20 flex gap-2">
               <input 
                 type="file" 
                 accept=".wav" 
                 ref={fileInputRef} 
                 className="hidden" 
                 onChange={handleFileUpload}
               />
               <button 
                 onClick={() => fileInputRef.current.click()} 
                 className="px-3 py-1.5 bg-slate-800 hover:bg-slate-700 border border-slate-700 rounded text-xs font-medium text-slate-300 flex items-center gap-2 shadow-lg"
               >
                 <Upload size={14}/> Open WAV
               </button>
               <button 
                 onClick={generateRandomTarget}
                 className="px-3 py-1.5 bg-slate-800 hover:bg-slate-700 border border-slate-700 rounded text-xs font-medium text-slate-300 flex items-center gap-2 shadow-lg"
               >
                 <Shuffle size={14}/> Random Target
               </button>
               <button 
                 onClick={runOptimization}
                 disabled={!measurement}
                 className={`px-3 py-1.5 rounded text-xs font-medium flex items-center gap-2 shadow-lg transition-colors
                   ${isOptimizing ? 'bg-rose-600 text-white animate-pulse' : 
                     !measurement ? 'bg-slate-800 text-slate-600 cursor-not-allowed border border-slate-800' : 'bg-cyan-600 text-white hover:bg-cyan-500'}
                 `}
               >
                 <Play size={14}/> {isOptimizing ? "Fitting..." : "Auto Fit"}
               </button>
            </div>
          )}

          <div className="flex-1 p-6 relative h-full flex flex-col">
            <div className="h-full w-full bg-slate-900/30 rounded-xl border border-slate-800/50 p-4 overflow-hidden flex flex-col">
              
              {viewMode === 'spl' && (
                <div className="flex flex-col h-full">
                  {/* Top: SPL Chart (40%) */}
                  <div className="h-[40%] min-h-0 relative border-b border-slate-800 pb-2">
                    <div className="absolute top-0 left-2 text-xs font-bold text-slate-500 z-10 pointer-events-none">Frequency Response (SPL)</div>
                    <ResponsiveContainer width="100%" height="100%">
                      <ComposedChart data={chartData} margin={{ top: 25, right: 30, bottom: 5, left: 10 }}>
                        <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#1e293b" />
                        <XAxis dataKey="frequency" type="number" scale="log" domain={[20, 20000]} ticks={[20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]} tickFormatter={(t) => t >= 1000 ? `${t/1000}k` : t} stroke="#475569" fontSize={10} tickMargin={5}/>
                        <YAxis domain={['auto', 'auto']} stroke="#475569" fontSize={10}>
                            <Label value="SPL (dB)" angle={-90} position="insideLeft" fill="#64748b" style={{ textAnchor: 'middle' }} />
                        </YAxis>
                        <Tooltip contentStyle={{ backgroundColor: '#0f172a', borderColor: '#334155', color: '#f1f5f9' }} labelFormatter={(l) => `${l} Hz`} formatter={(v) => v.toFixed(1) + ' dB'} />
                        <Legend verticalAlign="top" height={20} iconSize={10} wrapperStyle={{ fontSize: '10px' }}/>
                        {/* Target (Harman) - Aligned */}
                        <Line name="Target (Harman)" type="monotone" dataKey="target" stroke="#a21caf" strokeWidth={3} strokeDasharray="5 5" dot={false} activeDot={false} strokeOpacity={0.6} />
                        {/* Measurement (if any) */}
                        {measurement && (
                          <Line name="Measurement" type="monotone" dataKey="meas" stroke="#eab308" strokeWidth={2} strokeDasharray="3 3" dot={false} activeDot={false} />
                        )}
                        {/* Simulation */}
                        <Line name="Simulation" type="monotone" dataKey="spl" stroke="#06b6d4" strokeWidth={2} dot={false} activeDot={{ r: 4, fill: '#fff' }} />
                      </ComposedChart>
                    </ResponsiveContainer>
                  </div>

                  {/* Middle: IR Chart (35%) */}
                  <div className="h-[35%] min-h-0 relative pt-2 border-b border-slate-800 pb-2">
                    <div className="absolute top-2 left-2 text-xs font-bold text-slate-500 z-10 pointer-events-none">Impulse Response</div>
                    <ResponsiveContainer width="100%" height="100%">
                      <ComposedChart data={irData} margin={{ top: 25, right: 30, bottom: 20, left: 10 }}>
                        <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#1e293b" />
                        <XAxis dataKey="time" type="number" domain={[0, 3]} stroke="#475569" fontSize={10}>
                            <Label value="Time (ms)" position="bottom" offset={0} fontSize={10} fill="#64748b"/>
                        </XAxis>
                        <YAxis domain={[-1.1, 1.1]} stroke="#475569" fontSize={10}>
                            <Label value="Amp" angle={-90} position="insideLeft" fill="#64748b" style={{ textAnchor: 'middle' }} />
                        </YAxis>
                        <ReferenceLine x={0} stroke="#ef4444" strokeDasharray="3 3" />
                        <Scatter data={irPeak ? [irPeak] : []} fill="#ef4444" shape="circle" legendType="none" />
                        <Tooltip contentStyle={{ backgroundColor: '#0f172a', borderColor: '#334155', color: '#f1f5f9' }} labelFormatter={(l) => `${l} ms`} />
                        {measurement && (
                          <Line type="monotone" dataKey="measAmp" stroke="#eab308" strokeWidth={1.5} dot={false} strokeOpacity={0.7} />
                        )}
                        <Line type="monotone" dataKey="amplitude" stroke="#06b6d4" strokeWidth={1.5} dot={false} />
                      </ComposedChart>
                    </ResponsiveContainer>
                  </div>

                  {/* Bottom: Analysis Dashboard (Restored) */}
                  <div className="flex-1 flex pt-4 overflow-hidden">
                    {/* Tags */}
                    <div className="w-1/2 border-r border-slate-800 pr-4 flex flex-col">
                      <div className="flex items-center gap-2 mb-3">
                        <FileText className="w-4 h-4 text-slate-500" />
                        <h2 className="text-sm font-semibold text-slate-300">Timbre Analysis</h2>
                      </div>
                      <div className="flex flex-wrap gap-2 content-start overflow-y-auto">
                        {timbre.tags.length > 0 ? (
                          timbre.tags.map((tag, i) => (
                            <span key={i} className={`px-2 py-1 rounded-md text-xs font-medium border
                              ${tag.includes("Heavy") || tag.includes("Strong") || tag.includes("Boomy") ? "bg-rose-950/30 text-rose-400 border-rose-900/50" : 
                                tag.includes("Bright") || tag.includes("Airy") || tag.includes("Sharp") ? "bg-amber-950/30 text-amber-400 border-amber-900/50" :
                                tag.includes("Dark") || tag.includes("Recessed") || tag.includes("Lack") ? "bg-slate-800 text-slate-400 border-slate-700" :
                                "bg-cyan-950/30 text-cyan-400 border-cyan-900/50"
                              }`}>
                              {tag}
                            </span>
                          ))
                        ) : (
                          <span className="text-slate-500 text-xs italic">Neutral Timbre (Matches Harman Target)</span>
                        )}
                      </div>
                    </div>

                    {/* Spectrum Bars */}
                    <div className="w-1/2 pl-4 flex flex-col">
                      <div className="flex items-center gap-2 mb-3">
                        <BarChart2 className="w-4 h-4 text-slate-500" />
                        <h2 className="text-sm font-semibold text-slate-300">Spectrum Deviation</h2>
                      </div>
                      <div className="flex-1 flex flex-col justify-between overflow-y-auto">
                        {['subBass', 'bass', 'mid', 'treble', 'air'].map((band) => {
                          const val = timbre.analysis[band] || 0;
                          const absVal = Math.abs(val);
                          const barWidth = Math.min(absVal * 10, 100); 
                          const isPos = val > 0;
                          
                          return (
                            <div key={band} className="flex items-center gap-2 text-[10px]">
                              <span className="w-10 font-medium text-slate-500 uppercase text-right">{band}</span>
                              <div className="flex-1 h-3 bg-slate-800 rounded-full relative overflow-hidden flex items-center">
                                <div className="absolute left-1/2 w-px h-full bg-slate-600 z-10"></div>
                                <div 
                                  className={`h-full rounded-sm transition-all duration-300 ${Math.abs(val) > 2.5 ? (isPos ? 'bg-rose-500' : 'bg-indigo-500') : 'bg-slate-600'}`}
                                  style={{ 
                                    width: `${barWidth}%`, 
                                    marginLeft: isPos ? '50%' : `calc(50% - ${barWidth}%)`
                                  }}
                                ></div>
                              </div>
                              <span className={`w-8 font-mono text-right ${Math.abs(val) > 2.5 ? 'text-slate-200 font-bold' : 'text-slate-500'}`}>
                                {val > 0 ? '+' : ''}{val.toFixed(1)}
                              </span>
                            </div>
                          );
                        })}
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {viewMode === 'circuit' && <CircuitDiagram />}

            </div>
          </div>
        </div>
      </div>
    </div>
  );
}