/**
 * Bayes-Estimation: Clinical & Scientific Presets
 * Pre-configured real-world biomedical and neurosurgical data presets.
 */

export const Presets = {
  twoProportions: [
    {
      id: "stroke_trial",
      title: "Acute Ischemic Stroke: Endovascular Thrombectomy vs Medical Therapy",
      description: "Functional independence (mRS 0-2 at 90 days) comparing EVT vs best medical therapy.",
      kA: 58, nA: 267, // Control: 21.7%
      kB: 114, nB: 267, // Treatment: 42.7%
      priorA: 1, priorB: 1,
      ropeLow: -0.03, ropeHigh: 0.03
    },
    {
      id: "craniosynostosis",
      title: "Endoscopic vs Open Sagittal Synostosis: Blood Transfusion Rate",
      description: "Comparing intraoperative blood transfusion rates between minimally invasive endoscopic strip craniectomy vs open cranial vault remodeling.",
      kA: 42, nA: 50, // Open: 84.0%
      kB: 7, nB: 50, // Endoscopic: 14.0%
      priorA: 1, priorB: 1,
      ropeLow: -0.05, ropeHigh: 0.05
    },
    {
      id: "sah_vasospasm",
      title: "Aneurysmal SAH: Prophylactic Nimodipine vs Control",
      description: "Symptomatic delayed cerebral ischemia (vasospasm) in aneurysmal subarachnoid hemorrhage.",
      kA: 68, nA: 175, // Control: 38.9%
      kB: 38, nB: 179, // Nimodipine: 21.2%
      priorA: 1, priorB: 1,
      ropeLow: -0.04, ropeHigh: 0.04
    }
  ],

  twoMeans: [
    {
      id: "tbi_icp",
      title: "Severe TBI: Mean Intracranial Pressure (mmHg)",
      description: "Comparing intracranial pressure (ICP in mmHg) at 48 hours post-injury in Targeted Temperature Management vs Standard Normothermia.",
      meanA: 21.4, sdA: 4.8, nA: 64, // Standard
      meanB: 15.2, sdB: 3.9, nB: 62, // Targeted Temperature
      ropeLow: -1.0, ropeHigh: 1.0,
      priorScale: 0.707
    },
    {
      id: "inph_gait",
      title: "iNPH Tap Test: Timed-Up-and-Go (TUG) Seconds",
      description: "Time to complete 3-meter Timed-Up-and-Go test before vs after 40mL CSF drainage in normal pressure hydrocephalus.",
      meanA: 24.8, sdA: 5.6, nA: 38, // Baseline
      meanB: 18.2, sdB: 4.9, nB: 38, // Post-Tap
      ropeLow: -1.5, ropeHigh: 1.5,
      priorScale: 0.707
    },
    {
      id: "glioma_eor",
      title: "Glioma Resection: Extent of Resection (% Volume)",
      description: "Percentage volumetric tumor resection achieved with Intraoperative MRI vs Standard Neuronavigation.",
      meanA: 87.5, sdA: 8.4, nA: 52, // Standard
      meanB: 95.8, sdB: 4.1, nB: 48, // Intraop MRI
      ropeLow: -2.0, ropeHigh: 2.0,
      priorScale: 0.707
    }
  ],

  singleProportion: [
    {
      id: "microdiscectomy_recurrence",
      title: "Lumbar Microdiscectomy: 2-Year Recurrent Disc Herniation Rate",
      description: "Symptomatic re-herniation requiring revision surgery within 24 months.",
      k: 14, n: 210, // 6.7%
      priorA: 1, priorB: 15, // Skeptical historical prior ~ 6%
      benchmark: 0.10
    },
    {
      id: "shunt_infection",
      title: "Ventriculoperitoneal Shunt Infection Rate (Antibiotic-Impregnated Catheters)",
      description: "Shunt infection within 6 months of catheter implantation.",
      k: 5, n: 160, // 3.1%
      priorA: 1, priorB: 1,
      benchmark: 0.05
    }
  ],

  diagnosticNomogram: [
    {
      id: "inph_tap_test",
      title: "iNPH Diagnosis: CSF Tap Test Accuracy",
      description: "Clinical diagnostic accuracy of 40mL CSF Tap Test in predicting positive VP shunt response in suspected idiopathic NPH.",
      preTestProb: 0.40, // 40% clinical suspicion
      sensitivity: 0.65, // Tap test has moderate sensitivity
      specificity: 0.88, // Tap test has high specificity
      sampleN_disease: 85,
      sampleN_healthy: 70
    },
    {
      id: "carotid_stenosis",
      title: "Severe Carotid Stenosis (≥70%): Duplex Ultrasound PSV",
      description: "Peak Systolic Velocity > 230 cm/s for diagnosing severe internal carotid artery stenosis.",
      preTestProb: 0.25,
      sensitivity: 0.90,
      specificity: 0.92,
      sampleN_disease: 120,
      sampleN_healthy: 150
    }
  ]
};
