package opus

const OPUS_AUTO = -1000

// / <summary>
// / Maximum bitrate
// / </summary>
const OPUS_BITRATE_MAX = -1

// from analysis.c
const NB_FRAMES = 8
const NB_TBANDS = 18
const NB_TOT_BANDS = 21
const NB_TONAL_SKIP_BANDS = 9
const ANALYSIS_BUF_SIZE = 720

/* 15 ms at 48 kHz */
const DETECT_SIZE = 200

const MAX_ENCODER_BUFFER = 480
