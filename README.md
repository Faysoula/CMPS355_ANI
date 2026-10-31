# c4bot — local Connect 4 bot for papergames.io

Plays perfect Connect 4 against the papergames robot. Runs on your laptop, no
browser extension.

## Run it

1. Install Node.js (LTS) from https://nodejs.org if you don't have it.
2. Double-click `start.bat` (first run installs two packages, ~10 MB).
3. A Chrome window opens on papergames. Log in once if you want (the login is
   kept in `browser-profile/`), then click **Play vs robot**.
4. The bot plays every move for you. Watch the terminal for what it's doing.

Any color, either side starting — it works out whose turn it is from the board.

| Hotkey      | What it does                                        |
|-------------|-----------------------------------------------------|
| `Shift + 0` | pause / resume auto-play                            |
| `Shift + =` | play the best move right now (works while paused)   |
| `Shift + -` | show every column's score (+ win, 0 draw, − loss)   |

## How it works

- **`engine.mjs`** — the solver from
  [benjaminrall/connect-four-ai](https://github.com/benjaminrall/connect-four-ai):
  Rust compiled to WebAssembly (near-native speed) with an opening book.
  Every position is solved exactly, usually in under a millisecond, worst case
  a few seconds around moves 9–13. No time limits, no guessing.
- **`bot.mjs`** — opens Chrome (or Edge) with Playwright, reads the board
  from the page every 150 ms, waits until pieces stop falling, and clicks the
  best column with a real mouse click.
- Auto-play only runs in rooms with a robot player (`ROBOT_NAMES` at the top
  of `bot.mjs`, default `"Paper Man"`). The terminal prints the player names
  for every room, so if your robot level has a different name, add it there.
- Only the classic 7×6 board is supported.

## What to expect

Connect 4 is solved: whoever moves first wins with perfect play. So the bot
wins every game where it starts. When the robot starts, a perfect robot would
win, but the bot always picks the move that holds out longest, and any mistake
by the robot gets punished immediately.

## Test

`npm test` runs the real bot (headless) against `test/mock.html`, a local copy
of the papergames board markup with a robot opponent:

```
node test/run-mock.mjs 10 greedy     # weaker robot: bot should win all
node test/run-mock.mjs 10 perfect    # perfect robot: bot wins exactly when it starts
```

This needs a Chromium for Playwright: `npx playwright-core install chromium`.
