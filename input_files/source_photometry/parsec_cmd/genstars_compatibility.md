# genstars Compatibility Check

This compares the regenerated PARSEC/CMD normalized tables against the copied
`genstars` processed source-photometry tables. Differences are expected for the
low-mass empirical/Baraffe parts of the legacy tables, so the primary checks use
the PARSEC-dominated ranges.

### prime thin1

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 13 | 0.008 | 0.023 |
| MI_c vs Imag | 13 | 0.01 | 0.018 |
| MJ_2M vs Jmag_2mass | 13 | 0.007 | 0.018 |
| MH_2M vs Hmag_2mass | 13 | 0.032 | 0.039 |
| MK_2M vs Ksmag_2mass | 13 | 0.012 | 0.016 |

### roman thin1

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 29 | 0.004 | 0.016 |
| MH_2M vs Hmag_2mass | 29 | 0.003 | 0.014 |
| MK_2M vs Ksmag_2mass | 29 | 0.003 | 0.015 |
| MZ087 vs F087mag | 29 | 0.007 | 0.015 |
| MW146 vs F146mag | 29 | 0.003 | 0.015 |
| MF213 vs F213mag | 29 | 0.003 | 0.015 |

### prime thin2

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 12 | 0.008 | 0.009 |
| MI_c vs Imag | 12 | 0.002 | 0.003 |
| MJ_2M vs Jmag_2mass | 12 | 0.01 | 0.014 |
| MH_2M vs Hmag_2mass | 12 | 0.031 | 0.032 |
| MK_2M vs Ksmag_2mass | 12 | 0.01 | 0.013 |

### roman thin2

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 29 | 0 | 0 |
| MH_2M vs Hmag_2mass | 29 | 0 | 0 |
| MK_2M vs Ksmag_2mass | 29 | 0 | 0 |
| MZ087 vs F087mag | 29 | 0 | 0 |
| MW146 vs F146mag | 29 | 0 | 0 |
| MF213 vs F213mag | 29 | 0 | 0 |

### prime thin3

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 11 | 0.039 | 0.054 |
| MI_c vs Imag | 11 | 0.016 | 0.02 |
| MJ_2M vs Jmag_2mass | 11 | 0.019 | 0.027 |
| MH_2M vs Hmag_2mass | 11 | 0.02 | 0.028 |
| MK_2M vs Ksmag_2mass | 11 | 0.003 | 0.011 |

### roman thin3

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 29 | 0.011 | 0.067 |
| MH_2M vs Hmag_2mass | 29 | 0.009 | 0.069 |
| MK_2M vs Ksmag_2mass | 29 | 0.007 | 0.066 |
| MZ087 vs F087mag | 29 | 0.018 | 0.072 |
| MW146 vs F146mag | 29 | 0.01 | 0.067 |
| MF213 vs F213mag | 29 | 0.008 | 0.066 |

### prime thin4

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 10 | 0.06 | 0.079 |
| MI_c vs Imag | 10 | 0.027 | 0.031 |
| MJ_2M vs Jmag_2mass | 10 | 0.028 | 0.035 |
| MH_2M vs Hmag_2mass | 10 | 0.017 | 0.02 |
| MK_2M vs Ksmag_2mass | 10 | 0.003 | 0.01 |

### roman thin4

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 24 | 0.016 | 0.058 |
| MH_2M vs Hmag_2mass | 24 | 0.013 | 0.06 |
| MK_2M vs Ksmag_2mass | 24 | 0.011 | 0.056 |
| MZ087 vs F087mag | 24 | 0.028 | 0.066 |
| MW146 vs F146mag | 24 | 0.015 | 0.059 |
| MF213 vs F213mag | 24 | 0.012 | 0.056 |

### prime thin5

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 11 | 0.088 | 0.129 |
| MI_c vs Imag | 11 | 0.047 | 0.056 |
| MJ_2M vs Jmag_2mass | 11 | 0.039 | 0.051 |
| MH_2M vs Hmag_2mass | 11 | 0.008 | 0.011 |
| MK_2M vs Ksmag_2mass | 11 | 0.009 | 0.021 |

### roman thin5

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 25 | 0.028 | 0.047 |
| MH_2M vs Hmag_2mass | 25 | 0.021 | 0.049 |
| MK_2M vs Ksmag_2mass | 25 | 0.018 | 0.045 |
| MZ087 vs F087mag | 25 | 0.047 | 0.105 |
| MW146 vs F146mag | 25 | 0.027 | 0.056 |
| MF213 vs F213mag | 25 | 0.019 | 0.048 |

### prime thin6

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 13 | 0.184 | 0.227 |
| MI_c vs Imag | 13 | 0.106 | 0.125 |
| MJ_2M vs Jmag_2mass | 13 | 0.082 | 0.106 |
| MH_2M vs Hmag_2mass | 13 | 0.019 | 0.061 |
| MK_2M vs Ksmag_2mass | 13 | 0.032 | 0.078 |

### roman thin6

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 27 | 0.058 | 0.105 |
| MH_2M vs Hmag_2mass | 27 | 0.037 | 0.087 |
| MK_2M vs Ksmag_2mass | 27 | 0.031 | 0.083 |
| MZ087 vs F087mag | 27 | 0.095 | 0.15 |
| MW146 vs F146mag | 27 | 0.054 | 0.101 |
| MF213 vs F213mag | 27 | 0.032 | 0.083 |

### prime thin7

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 12 | 0.265 | 0.321 |
| MI_c vs Imag | 12 | 0.146 | 0.283 |
| MJ_2M vs Jmag_2mass | 12 | 0.112 | 0.274 |
| MH_2M vs Hmag_2mass | 12 | 0.049 | 0.237 |
| MK_2M vs Ksmag_2mass | 12 | 0.063 | 0.255 |

### roman thin7

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 24 | 0.083 | 0.273 |
| MH_2M vs Hmag_2mass | 24 | 0.059 | 0.263 |
| MK_2M vs Ksmag_2mass | 24 | 0.05 | 0.26 |
| MZ087 vs F087mag | 24 | 0.135 | 0.284 |
| MW146 vs F146mag | 24 | 0.074 | 0.271 |
| MF213 vs F213mag | 24 | 0.052 | 0.26 |

### prime bar

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 30 | 0.008 | 0.009 |
| MI_c vs Imag | 30 | 0.003 | 0.003 |
| MJ_2M vs Jmag_2mass | 30 | 0.004 | 0.015 |
| MH_2M vs Hmag_2mass | 30 | 0.028 | 0.031 |
| MK_2M vs Ksmag_2mass | 30 | 0.008 | 0.013 |

### roman bar

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 47 | 0 | 0 |
| MH_2M vs Hmag_2mass | 47 | 0 | 0 |
| MK_2M vs Ksmag_2mass | 47 | 0 | 0 |
| MZ087 vs F087mag | 47 | 0 | 0 |
| MW146 vs F146mag | 47 | 0 | 0 |
| MF213 vs F213mag | 47 | 0 | 0 |

### prime NSD

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 12 | 0.184 | 0.248 |
| MI_c vs Imag | 12 | 0.119 | 0.135 |
| MJ_2M vs Jmag_2mass | 12 | 0.063 | 0.094 |
| MH_2M vs Hmag_2mass | 12 | 0.07 | 0.099 |
| MK_2M vs Ksmag_2mass | 12 | 0.042 | 0.071 |

### roman NSD

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 26 | 0.055 | 0.11 |
| MH_2M vs Hmag_2mass | 26 | 0.049 | 0.09 |
| MK_2M vs Ksmag_2mass | 26 | 0.042 | 0.084 |
| MZ087 vs F087mag | 26 | 0.106 | 0.153 |
| MW146 vs F146mag | 26 | 0.06 | 0.104 |
| MF213 vs F213mag | 26 | 0.044 | 0.084 |

### prime thick

Mass range: `0.55 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MV_j vs Vmag | 73 | 0.156 | 1.267 |
| MI_c vs Imag | 73 | 0.183 | 1.376 |
| MJ_2M vs Jmag_2mass | 73 | 0.218 | 1.453 |
| MH_2M vs Hmag_2mass | 73 | 0.268 | 1.565 |
| MK_2M vs Ksmag_2mass | 73 | 0.251 | 1.553 |

### roman thick

Mass range: `0.09 <= Mini <= 1.0` and `abs(MPD - Mini) <= 0.01`.

| column | n | median_abs_diff | max_abs_diff |
| --- | ---: | ---: | ---: |
| MJ_2M vs Jmag_2mass | 87 | 0.116 | 1.459 |
| MH_2M vs Hmag_2mass | 87 | 0.125 | 1.537 |
| MK_2M vs Ksmag_2mass | 87 | 0.127 | 1.546 |
| MZ087 vs F087mag | 87 | 0.107 | 1.389 |
| MW146 vs F146mag | 87 | 0.119 | 1.486 |
| MF213 vs F213mag | 87 | 0.127 | 1.546 |

