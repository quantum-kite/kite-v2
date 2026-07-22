# KITE Documentation — Ultimate Goal

This is the north star for every stage of this project (discussion, execution, audit).
Check new work against this before considering it done — not just "is this accurate,"
but "does this actually move a real user closer to this outcome."

## The goal

A scientist who has never used KITE before should be able to:

1. **Download and install it with minimal friction.** No guessing dependencies, no
   figuring out undocumented build flags, no trial and error. If something about
   installation is currently fragile or under-explained, that's as much a priority
   as any theory or API gap.

2. **Actually understand what they're running, not just how to invoke it.** Usage
   documentation should explain the physics behind each calculation (what a
   quantity means, why the method works, what approximations are involved) *and*
   the functions/parameters well enough that a user can connect "this knob" to
   "this physical effect" — not just copy a working example without understanding it.

Both halves matter together. Easy installation with opaque usage just gets someone
running a black box. Deep physics explanation with a painful install means nobody
gets far enough to read it.

## How to use this file

- **In Plan/discussion:** when scoping any doc or install-process change, ask whether
  it's actually closing a gap in this goal, not just adding more content.
- **In execution (careful-executor):** when writing a physics explanation, connect it
  explicitly to the function/parameter that exercises it, and vice versa — don't leave
  physics and API documentation as disconnected sections a reader has to cross-reference
  themselves.
- **In audit (structure-auditor):** when reviewing a completed change, check it against
  this goal directly, not only against internal consistency — a technically correct,
  well-organized page that still doesn't help a new user install or understand KITE
  hasn't actually succeeded.
