# Rate control prototype / first version

## Simulation parameters

- per simulation:
  - `max_absolute_rate` - optional, `absolute_rate` gets capped at this value, in units of reactions per D/M step
  - `smoothing_constant` - mandatory with default values, two small numbers (between 0 and 1, close to 0), the smaller both are the less sensitive the smoothing is to noise, but the slower it will react to change
  - `highest_probability` - mandatory with default value 1 currently, a number between 0 and 1, the probability of a reaction happening for every reaction will get multiplied by this number, so it can be used to scale all reactions.
- per reaction:
  - `relative_rate` - rate constant `k` for the reaction, relative to `absolute_rate`, therefore contains no units of time, for bimolecular or larger order reactions: concentration units of molecules per simulation box

## State kept

- per simulation:
  - `absolute_rate` - multiply this by each `relative_rate` to get the `k` rate constant that the algorithm aims for
- per reaction:
  - `observed_rate` - connected to the frequency of passing geometry conditions of a reaction

## Goals and limitations

- Non-goal: realistic absolute rates
- Goal: correct instantaneous relative rates between the reactions in the system at a specific point of the simulation
- Built on the assumption, that the system contains all reactants evenly mixed, and that all matches of a single class of reactant (a single graph) react at roughly similar rates
- Limitations: `absolute_rate` can differ during the reaction, so the rate of a reaction can change (even by orders of magnitude) during the same simulation, for example if new, much much slower reactions start happening during the simulation. Generally this rate limiting is most useful if the nature of the system (what types of reactions happen) doesn't change that much during a simulation from the initial point. Or in some cases, it can mean that the initial rates will be different from most of the simulation (e.g. dimerization and breakage - before there is anything to break, we can't correctly set the absolute rate of dimerization to get the good relative rate with breakage, if breakage is much slower).
- Limitation: relative rates are meaningless with only a single reaction. Beware of unimolecular reactions with only a rate condition - if there is only a single breaking reaction in the system it will proceed at 100% speed (all done in a single D/M run) if `max_absolute_rate` or `highest_probability` are not set to limit it.
- Note: Uni and bimolecular reaction actual relative rates to eachother will be concentration dependent, since a dimensional analysis will reveal that `absolute_rate` has a unit problem
  - but this is expected, since saying that reaction A (unimolecular) should be 50x faster than reaction B (bimolecular) at all reactant concentrations is impossible
- currently all reactions in a simulation get tied to the same `absolute_rate`, it's not possible to "group them" (e.g. reactions A and B have a specific relative rate to eachother, independent from reactions C and D), this could be lifted later, maybe?

## Algorithm

- T* contains a value `absolute_rate`, which gives the dimension
  of rate constants (to get a real `k` for every reaction)
- initially (first D/M run) it is not set, the first D/M run is considered a warm up and no rate controlled reactions will happen
- After each D/M run, `observed_rate` is calculated:
  - as reactions during the last step, divided by `relative_rate` and the product of all concentrations (in units of molecules per simulation box) - this way it is converted into the units of `absolute_rate` and already scaled by the `relative_rate`!
  - the `observed_rate` value is kept smoothed using a double exponential smoothing (very simple smoothing formula that smoothes the values and slopes for values)
  - if there are no reactants for the reaction, its `observed_rate` is not updated
  - initially `observed_rate` is not set
- `absolute_rate` is set to the smallest of `observed_rate` * `highest_probability`, or `max_absolute_rate` (if specified, and if it is smaller than all of them)
- in the detection algorithm, if a reaction is rate controlled but `absolute_rate` or the reaction specific `observed_rate` are not set, the reaction is rejected
- in the detection algorithm, the probability of accepting the reaction is `absolute_rate` / `observed_rate`
  - Note: for the slowest reaction, `absolute_rate` = `observed_rate`, therefore this probability will be 100% (unless `max_absolute_rate` is smaller)
  - for every other reaction, `observed_rate` is bigger than `absolute_rate`, by a factor that estimates the relative frequency of geometry conditions becoming true between the two reactions
