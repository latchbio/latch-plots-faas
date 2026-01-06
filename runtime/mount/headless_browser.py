"""Minimal Playwright wrapper for headless notebook browser."""

from playwright.async_api import async_playwright, Browser, Page, Playwright


class HeadlessBrowser:
    """Manages a headless Chromium browser for running notebook frontend."""

    def __init__(self) -> None:
        self.playwright: Playwright | None = None
        self.browser: Browser | None = None
        self.page: Page | None = None

    async def start(self, notebook_url: str) -> Page:
        """Launch browser and navigate to notebook.

        Args:
            notebook_url: Full URL to the notebook page with auth params.

        Returns:
            The Page object for the loaded notebook.
        """
        print(f"[headless] Launching Chromium...")
        self.playwright = await async_playwright().start()
        self.browser = await self.playwright.chromium.launch(headless=True)
        self.page = await self.browser.new_page(viewport={"width": 1280, "height": 800})
        print(f"[headless] Navigating to {notebook_url[:100]}...")
        await self.page.goto(notebook_url, wait_until="networkidle")
        print(f"[headless] Page loaded, URL: {self.page.url[:100]}")
        return self.page

    async def wait_for_agent_ready(self, timeout_ms: int = 30000) -> None:
        """Wait for frontend to signal agent is connected.

        Args:
            timeout_ms: Maximum time to wait in milliseconds.

        Raises:
            TimeoutError: If agent doesn't become ready within timeout.
        """
        if self.page is None:
            raise RuntimeError("Browser not started. Call start() first.")
        print(f"[headless] Waiting for [data-agent-ready='true'] (timeout={timeout_ms}ms)...")
        try:
            await self.page.wait_for_selector("[data-agent-ready='true']", timeout=timeout_ms)
            print("[headless] Agent ready!")
        except Exception as e:
            # Log current page state for debugging
            try:
                agent_ready_attr = await self.page.evaluate(
                    "document.querySelector('[data-agent-ready]')?.getAttribute('data-agent-ready')"
                )
                print(f"[headless] Timeout waiting for agent ready. Current data-agent-ready={agent_ready_attr}")
            except:
                print(f"[headless] Timeout waiting for agent ready. Could not read attribute.")
            raise

    async def stop(self) -> None:
        """Close browser and clean up resources."""
        if self.browser:
            await self.browser.close()
            self.browser = None
        if self.playwright:
            await self.playwright.stop()
            self.playwright = None
        self.page = None

