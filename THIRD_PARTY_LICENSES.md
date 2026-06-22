# Third-Party Components

CADENCE bundles a small number of third-party helper functions, listed below.
Each remains under its original license; CADENCE’s own MIT license does not
override them. All are permissive (BSD-2-Clause or public domain) and compatible
with CADENCE’s MIT distribution, provided each component’s copyright/attribution
notice is retained in its source file (it is).

| Component | File | Author | Source | License | Notice status |
|---|---|---|---|---|---|
| contourcs | `arrhythmia_dynamics_helper/contourcs.m` | Takeshi Ikuma (2010) | File Exchange | BSD-2-Clause | ✅ full text in file |
| real2rgb | `utils/real2rgb.m` | Oliver Woodford (2009–2010) | File Exchange | BSD-2-Clause | ✅ copyright + source/license pointer |
| rescale (→ `rescale_sat`) | `utils/rescale_sat.m` | Oliver Woodford (2009–2011) | File Exchange | BSD-2-Clause | ✅ copyright + pointer; renamed & lightly modified for CADENCE |
| rgb | `utils/rgb.m` | Kristján Jónasson, Univ. of Iceland (2009) | [File Exchange #24497](https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2) | Public domain | ✅ complete in file |
| DiscreteFrechetDist | `arrhythmia_dynamics_helper/frechet/DiscreteFrechetDist.m` | Zachary Danziger (2011) | [File Exchange #31922](https://www.mathworks.com/matlabcentral/fileexchange/31922-discrete-frechet-distance) | BSD-2-Clause | ✅ attribution restored (was initials only) |

## Licenses

### Public domain — `rgb`
`rgb.m` is released into the public domain by its author (Kristján Jónasson) and
carries no conditions. Attribution is retained in the file as a courtesy.

### BSD 2-Clause — `contourcs`, `real2rgb`, `rescale_sat`, `DiscreteFrechetDist`
These MATLAB File Exchange functions are provided under the BSD 2-Clause license.
Copyright is held by their respective authors (Takeshi Ikuma; Oliver Woodford;
Zachary Danziger), as noted in each file.

```
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Notes

- **Final license confirmation.** `contourcs` carries the full BSD text inline, so
  it is confirmed. For `real2rgb` / `rescale` (Woodford) and `DiscreteFrechetDist`
  (Danziger) — all older submissions — the BSD-2-Clause designation follows the
  MATLAB File Exchange standard; do a final eyeball of each submission’s
  “View License” link before public release. `rgb` states its public-domain status
  in-file, so it is confirmed.

- **MATLAB Runtime** — the compiled CADENCE application requires the MATLAB
  Runtime (© The MathWorks, Inc.), redistributed under MathWorks’ terms. It is a
  separately installed runtime and is not part of CADENCE.

- **License compatibility** — BSD-2-Clause, BSD-3-Clause, MIT, and public-domain
  components all combine freely with CADENCE’s MIT license.

- **Minor cleanup** — `rgb.m`’s author name contains mis-encoded characters
  (should read “Kristján Jónasson”); re-save the file as UTF-8 to fix. Cosmetic,
  not a compliance issue.

_Last reviewed: 2026-06-21._
