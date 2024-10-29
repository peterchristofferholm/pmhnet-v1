DROP SCHEMA IF EXISTS pmhnet CASCADE;

CREATE SCHEMA pmhnet
    CREATE TABLE surv_master (
        "pid"        INTEGER PRIMARY KEY,
        "train"      BOOLEAN,     -- split indicator
        "time"       INTEGER,     -- survival time in days
        "event"      BOOLEAN      -- 1 = death, 0 = censoring
    )
    CREATE TABLE baseline (
        "pid"        INTEGER UNIQUE REFERENCES surv_master (pid),
        "sex"        BOOLEAN,     -- 1 = male, 0 = female
        "age"        FLOAT8,      -- age at debut
        "height"     NUMERIC(5),  -- height at baseline (cm)
        "weight"     NUMERIC(5),  -- weight at baseline (kg)
        "pulse"      INTEGER,     -- heart rate (bpm)
        "sys_bp"     INTEGER,     -- systolic blood pressure (mmHg)
        "dia_bp"     INTEGER,     -- diastolic blood pressure (mmHg)
        "crea"       NUMERIC(5),  -- serum creatinine (umol/L)
        "arrest"     BOOLEAN,     -- cardiac arrest at admission (0/1)
        "enzymes"    BOOLEAN,     -- increased cardiac enzymes (0/1)
        "stemi"      BOOLEAN,     -- stemi (0/1)
        "smoking"    VARCHAR(3)
            CHECK (smoking IN ('yes', 'ex', 'no')),
        "vessels"    VARCHAR(3)   -- number of affected vessels
            CHECK (vessels IN ('DIF', '1VD', '2VD', '3VD')),
        "dominance"  VARCHAR(1)
            CHECK (dominance IN ('R', 'L', 'B')),
        "nyha"       INTEGER
            CHECK (nyha BETWEEN 1 AND 4),
        "ccs"        INTEGER
            CHECK (ccs BETWEEN 0 AND 4),
        "killip"     INTEGER
            CHECK (killip BETWEEN 1 and 4)
    )
    CREATE TABLE diagnoses (
        "pid"        INTEGER REFERENCES surv_master (pid),
        "time"       INTEGER,  -- time of diagnosis relative to index
        "level"      VARCHAR(5)
            CHECK ("level" in ('lvl4', 'lvl3', 'block', 'chap')),
        "code"       VARCHAR(8)
    )
    CREATE TABLE biochem (
        "pid"        INTEGER REFERENCES surv_master (pid),
        "qid"        VARCHAR NOT NULL,  -- quantity identifier
        "time"       INTEGER NOT NULL,  -- time of test relative to index
        "value"      INTEGER NOT NULL  -- "flagging" system (-1/0/1)
    )
    CREATE TABLE risk_scores (
        "pid"        INTEGER REFERENCES surv_master (pid),
        "time"       INTEGER,  -- prediction time (days)
        "type"       VARCHAR,
        "imputed"    BOOLEAN,  -- any grace values imputed? 
        "value"      FLOAT8
    )
    CREATE INDEX pid_idx_diag ON diagnoses (pid)
    CREATE INDEX pid_idx_base ON baseline (pid)
;
