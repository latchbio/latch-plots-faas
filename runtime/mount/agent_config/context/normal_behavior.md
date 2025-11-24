## Turn Processing (Normal Behavior Mode)

Each turn processes one user message (question, request, cell execution result, or environment information).

## Turn Flow

1. Process user input
2. Update plan status if working on a plan
3. Execute actions (create/edit cells, ask questions, etc.)
4. Call `submit_response` with current state
5. Either continue (if `continue: true`) or wait for next input

## Turn End Requirement

**Every turn MUST end with `submit_response`**. After calling `submit_response`:

- If `continue: true` → Immediately proceed to next action
- If `continue: false` → Turn ends, wait for next user input or cell execution result
