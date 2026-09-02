# Search strategies (verbatim)

Searches executed on 29 July 2026, covering 1995 to 2026, restricted to English.

## Bibliometric component

### Web of Science Core Collection
Field tag TS (topic). Refined to document types Article and Review.
```
TS=(anthocyanin* OR anthocyanidin* OR cyanidin OR delphinidin OR malvidin OR
pelargonidin OR peonidin OR petunidin OR "anthocyanin-rich")
AND
TS=(hypertension OR "blood pressure" OR "arterial pressure" OR antihypertensive
OR "vascular function" OR "endothelial function" OR "arterial stiffness" OR
"pulse wave velocity" OR "flow-mediated dilation")
```

### Scopus
```
TITLE-ABS-KEY(anthocyanin* OR anthocyanidin* OR cyanidin OR delphinidin OR
malvidin OR pelargonidin OR peonidin OR petunidin OR "anthocyanin-rich")
AND TITLE-ABS-KEY(hypertension OR "blood pressure" OR "arterial pressure" OR
antihypertensive OR "vascular function" OR "endothelial function" OR
"arterial stiffness" OR "pulse wave velocity" OR "flow-mediated dilation")
AND (LIMIT-TO(DOCTYPE,"ar") OR LIMIT-TO(DOCTYPE,"re"))
AND (LIMIT-TO(LANGUAGE,"English"))
AND PUBYEAR > 1994 AND PUBYEAR < 2027
```

### PubMed
```
(anthocyanin*[tiab] OR anthocyanidin*[tiab] OR cyanidin[tiab] OR
delphinidin[tiab] OR malvidin[tiab] OR pelargonidin[tiab] OR peonidin[tiab] OR
petunidin[tiab] OR "Anthocyanins"[Mesh])
AND
(hypertension[tiab] OR "blood pressure"[tiab] OR "arterial pressure"[tiab] OR
antihypertensive[tiab] OR "vascular function"[tiab] OR "endothelial
function"[tiab] OR "arterial stiffness"[tiab] OR "pulse wave velocity"[tiab] OR
"flow-mediated dilation"[tiab] OR "Hypertension"[Mesh] OR "Blood Pressure"[Mesh])
```

## Meta-analytic component

### PubMed, with the Cochrane randomized trial filter
```
[population block as above]
AND [outcome block as above]
AND (randomized controlled trial[pt] OR controlled clinical trial[pt] OR
randomized[tiab] OR randomised[tiab] OR placebo[tiab] OR "clinical trials as
topic"[Mesh] OR randomly[tiab] OR trial[ti])
NOT (animals[mh] NOT humans[mh])
```

### Cochrane CENTRAL
Trials only; no design filter required.
```
#1 (anthocyanin* OR anthocyanidin* OR cyanidin OR delphinidin OR malvidin OR
   pelargonidin OR peonidin OR petunidin):ti,ab,kw
#2 (hypertension OR "blood pressure" OR "arterial pressure" OR antihypertensive
   OR "vascular function" OR "endothelial function" OR "arterial stiffness" OR
   "pulse wave velocity" OR "flow-mediated dilation"):ti,ab,kw
#3 #1 AND #2, limited to Trials, 1995-2026
```

Embase was not accessible to the review team. Reference lists of included trials
and of relevant reviews were screened; ClinicalTrials.gov and the WHO ICTRP were
consulted for trials without a locatable publication.
