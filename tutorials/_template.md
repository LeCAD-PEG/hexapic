# A title that says what the reader will have done

<!-- Copyright 2026 David Lajevec — HEXAPIC
     The template of a tutorial. conf.py excludes it from the build. Copy it
     to tutorials/<name>.md and add <name> to the toctree in
     tutorials/index.md. -->

<div class="hx-meta"><span><b>Time</b> 30 min</span><span><b>Needs</b> a GPU, HEXAPIC built</span><span><b>Level</b> first steps</span></div>

:::{admonition} You will
:class: tip

- write a case for …
- run it and follow it in HEXAPIC / VIZ
- measure … and compare it with theory
:::

## 1. Set up the case

Explain what the case is and why it is interesting. Then give the commands
exactly as they run.

```bash
python python/examples/<name>.py runs/<name>
build/HEXAPIC runs/<name>/case.toml --check
```

## 2. Run it

```bash
build/HEXAPIC runs/<name>/case.toml
```

Say what to look at and what it should look like.

### Exercise

Change … and run again. What happens to …, and why?

:::{dropdown} Solution
Give the answer with the numbers that a run produces.
:::

## 3. Look at the result

## Summary

:::{admonition} Key points
:class: note

- …
:::

Next: {doc}`index`.
